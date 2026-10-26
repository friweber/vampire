//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>
#include <sstream>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "integrate.hpp"
#include "material.hpp"
#include "program.hpp"
#include "quantum.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      // reduced temperature per Kelvin for the spin system, kB/(hbar gamma_e)
      static double spin_temperature_scale(){
         return constants::kB/(constants::hbar*constants::gamma_SI);
      }

      // hbar in eV ps, the unit system of the lattice integrator
      static double hbar_eVps(){
         return constants::hbar*constants::kB_eV/constants::kB*1.0e12;
      }

      //---------------------------------------------------------------------------
      // Function to return the number of atoms this rank integrates
      //---------------------------------------------------------------------------
      int num_local_atoms(){
         #ifdef MPICF
            return vmpi::num_core_atoms + vmpi::num_bdry_atoms;
         #else
            return atoms::num_atoms;
         #endif
      }

      //---------------------------------------------------------------------------
      // Function to return the number of integration steps of the selected
      // program, or false if the program is not known to the module
      //---------------------------------------------------------------------------
      bool run_length(uint64_t& n_steps){

         const uint64_t et = sim::equilibration_time;
         const uint64_t tt = sim::total_time;
         const uint64_t lt = sim::loop_time;

         using program::program_t;
         switch(program::program){
            case program::program_t::benchmark: n_steps = tt;      return true; // benchmark
            case program::program_t::time_series: n_steps = et + tt; return true; // time series
            case program::program_t::hysteresis_loop: n_steps = et + lt; return true; // hysteresis
            case program::program_t::static_hysteresis_loop: n_steps = et + lt; return true; // static hysteresis
            case program::program_t::curie_temperature: n_steps = et + lt; return true; // Curie temperature
            case program::program_t::field_cool: n_steps = et + tt; return true; // field cool
            case program::program_t::laser_pulse: n_steps = et + tt; return true; // laser / temperature pulse
            case program::program_t::hamr_simulation: n_steps = et + tt; return true; // HAMR
            case program::program_t::lagrange_multiplier: n_steps = et + tt; return true; // Lagrange multiplier
            case program::program_t::partial_hysteresis_loop: n_steps = et + lt; return true; // partial hysteresis
            case program::program_t::localised_temperature_pulse: n_steps = et + tt; return true; // localised temperature pulse
            case program::program_t::effective_damping: n_steps = et + tt; return true; // effective damping
            case program::program_t::fmr: n_steps = et + tt; return true; // FMR
            case program::program_t::localised_field_cool: n_steps = et + tt; return true; // local field cool
            case program::program_t::electrical_pulse: n_steps = et + tt; return true; // electrical pulse
            case program::program_t::field_pulse: n_steps = et + tt; return true; // field pulse
            case program::program_t::domain_wall: n_steps = et + tt; return true; // domain walls
            case program::program_t::field_sweep: n_steps = et + lt; return true; // field sweep
            case program::program_t::spin_waves: n_steps = et + tt; return true; // spin waves
            default: return false;
         }

      }

      //---------------------------------------------------------------------------
      // Function to report whether the program moves sim::temperature during
      // the run; pre-generated noise is shaped once and cannot follow it
      //---------------------------------------------------------------------------
      bool dynamic_temperature_program(){

         using program::program_t;
         switch(program::program){
            case program::program_t::field_cool: // field cool
            case program::program_t::laser_pulse: // laser / temperature pulse
            case program::program_t::hamr_simulation: // HAMR
            case program::program_t::localised_temperature_pulse: // localised temperature pulse
            case program::program_t::localised_field_cool: // local field cool
               return true;
            default:
               return false;
         }

      }

      //---------------------------------------------------------------------------
      // Names for the log and the export header
      //---------------------------------------------------------------------------
      const char* spectrum_name(const spectrum_t s){
         if(s == classical) return "classical";
         if(s == quantum_zero) return "quantum";
         return "quantum-no-zero";
      }

      const char* generator_name(const generator_t g){
         return (g == on_the_fly) ? "on-the-fly" : "pre-generated";
      }

      //---------------------------------------------------------------------------
      // Function to check that a pre-generated bath can serve this program and
      // to return the run length used to size its window
      //---------------------------------------------------------------------------
      static uint64_t pre_generated_run_length(){

         uint64_t n_run = 0;

         if(noise_generator != pre_generated) return n_run;

         if(!run_length(n_run)){
            err::fatal("Error - program " + std::to_string(static_cast<int>(program::program)) +
                       " is not supported with pre-generated quantum noise. Use quantum:noise-generator = on-the-fly.");
         }

         if(dynamic_temperature_program()){
            err::fatal("Error - pre-generated quantum noise is shaped once, at the starting temperature, "
                       "and cannot follow a program that changes the temperature. "
                       "Use quantum:noise-generator = on-the-fly.");
         }

         return n_run;

      }

      //---------------------------------------------------------------------------
      // Thermostat (llg-quantum): per-material Lorentzian response and a spin
      // bath whose amplitude balances that response. Classical noise is
      // allowed here; it is the classical open-system model.
      //
      // Combined thermostat: with material:quantum-ohmic-damping alpha_e > 0
      // the spin is coupled to a second, Ohmic bath as well (Gilbert damping
      // alpha_e and its noise field, see ho_ohmic_equation_of_motion), with the
      // same spectrum and temperature but independent samples. Read as a
      // Lorentzian phonon bath plus an Ohmic electron bath. Its amplitude
      // sqrt(alpha_e/S0) is the one llg-heun uses for a coloured bath.
      //---------------------------------------------------------------------------
      static void initialize_thermostat(){

         const int n_mat = mp::num_materials;
         if(mp.size() < static_cast<size_t>(n_mat)) mp.resize(n_mat);

         material_A.resize(n_mat);
         material_gamma.resize(n_mat);
         material_omega0.resize(n_mat);
         material_S0.resize(n_mat);
         material_alpha_e.assign(n_mat, 0.0);

         std::vector<double> amp(n_mat);
         std::vector<double> coupling(n_mat);
         std::vector<double> amp_e(n_mat, 0.0);
         ohmic = false;

         for(int m = 0; m < n_mat; m++){

            const double alpha = mp::material[m].alpha;
            const double gamma = mp[m].gamma.get();
            const double omega0 = mp[m].omega0.get();

            if(gamma <= 0.0 || omega0 <= 0.0){
               err::fatal("Error - material " + std::to_string(m + 1) + " needs quantum-lorentzian-width and "
                          "quantum-lorentzian-central-frequency > 0 for sim:integrator = llg-quantum.");
            }

            // the moment in units of hbar*gamma_e converts the reduced temperature to a field
            const double S0 = mp::material[m].mu_s_SI/(constants::hbar*constants::gamma_SI);
            const double A = alpha*std::pow(omega0, 4)/gamma;

            material_A[m] = A;
            material_gamma[m] = gamma;
            material_omega0[m] = omega0;
            material_S0[m] = S0;

            // on-the-fly: the oscillator filters the bath, amplitude sqrt(gamma A/S0);
            // pre-generated: the filter is in the spectrum, amplitude 1/sqrt(S0)
            amp[m] = (noise_generator == pre_generated) ? 1.0/std::sqrt(S0) : std::sqrt(gamma*A/S0);
            coupling[m] = gamma*A;

            // Ohmic bath of the combined thermostat
            const double alpha_e = mp[m].alpha_e.get();
            material_alpha_e[m] = alpha_e;
            amp_e[m] = std::sqrt(alpha_e/S0);
            if(alpha_e > 0.0) ohmic = true;

         }

         if(ohmic && noise_generator == pre_generated){
            err::fatal("Error - material:quantum-ohmic-damping (combined thermostat) is supported with "
                       "quantum:noise-generator = on-the-fly only.");
         }

         const uint64_t n_run = pre_generated_run_length();
         const int n = num_local_atoms();

         setup_bath(spin_bath, noise_type, noise_generator, n, mp::dt, spin_temperature_scale(), amp, amp, true, n_run);
         spin_bath.coupling = coupling;
         if(ohmic){
            setup_bath(ohmic_bath, noise_type, on_the_fly, n, mp::dt, spin_temperature_scale(), amp_e, amp_e, false, 0);
            ohmic_bath.coupling = material_alpha_e;
         }
         allocate_thermostat(n);

         mode = thermostat;

         return;

      }

      //---------------------------------------------------------------------------
      // Spin bath for llg-heun: the classical thermal-field prefactor H_th_sigma
      // times sqrt(dt hbar gamma_e/2kB) reproduces H_th_sigma sqrt(T) in the
      // white limit
      //---------------------------------------------------------------------------
      static void initialize_llg_bath(){

         const int n_mat = mp::num_materials;
         const double pref = std::sqrt(mp::dt*constants::hbar*constants::gamma_SI/(2.0*constants::kB));

         std::vector<double> amp(n_mat);
         for(int m = 0; m < n_mat; m++) amp[m] = mp::material[m].H_th_sigma*pref;

         const uint64_t n_run = pre_generated_run_length();

         setup_bath(spin_bath, noise_type, noise_generator, num_local_atoms(), mp::dt, spin_temperature_scale(), amp, amp, false, n_run);

         mode = bath;

         return;

      }

      //---------------------------------------------------------------------------
      // Spin and lattice baths for spin-lattice dynamics. The lattice bath
      // runs in the lattice integrator's units (ps, eV) and reproduces the
      // classical force noise sqrt(2 eta kB T/m) in the white limit; both
      // baths switch to the equilibration amplitudes with the integrator.
      //---------------------------------------------------------------------------
      static void initialize_lattice_baths(){

         if(!lattice_parameters_set){
            err::fatal("Programmer error - quantum noise for spin-lattice dynamics requires the spin-lattice module "
                       "to be initialised first.");
         }

         const int n_mat = mp::num_materials;
         const int n = num_local_atoms();

         // spin bath
         const double pref = std::sqrt(mp::dt*constants::hbar*constants::gamma_SI/(2.0*constants::kB));
         std::vector<double> spin_amp(n_mat);
         std::vector<double> spin_amp_eq(n_mat);
         for(int m = 0; m < n_mat; m++){
            spin_amp[m] = mp::material[m].H_th_sigma*pref;
            spin_amp_eq[m] = mp::material[m].H_th_sigma_eq*pref;
         }

         const uint64_t n_run = pre_generated_run_length();

         setup_bath(spin_bath, noise_type, noise_generator, n, mp::dt, spin_temperature_scale(), spin_amp, spin_amp_eq, false, n_run);

         // lattice bath, always on-the-fly
         const double hbar_ps = hbar_eVps();
         std::vector<double> lattice_amp(n_mat);
         std::vector<double> lattice_amp_eq(n_mat);
         for(int m = 0; m < n_mat; m++){
            const double mass = lattice_mass[m];
            lattice_amp[m] = (mass > 0.0) ? std::sqrt(lattice_damping[m]*hbar_ps/mass) : 0.0;
            lattice_amp_eq[m] = (mass > 0.0) ? std::sqrt(lattice_damping_eq[m]*hbar_ps/mass) : 0.0;
         }

         setup_bath(lattice_bath, noise_type, on_the_fly, n, mp::dt_SI*1.0e12, constants::kB_eV/hbar_ps, lattice_amp, lattice_amp_eq, false, 0);

         mode = bath;

         return;

      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Function to initialise the quantum module
   //---------------------------------------------------------------------------
   void initialize(){

      using namespace internal;
      using integrate::integrator_t;

      mode = inactive;

      const bool coloured = (noise_type != classical);

      switch(integrate::integrator){

         case integrator_t::llg_quantum:
            initialize_thermostat();
            break;

         case integrator_t::llg_heun:
            if(coloured) initialize_llg_bath();
            break;

         case integrator_t::suzuki_trotter:
            if(coloured) initialize_lattice_baths();
            break;

         // Monte Carlo has no noise to colour
         case integrator_t::monte_carlo:
         case integrator_t::cmc:
         case integrator_t::hybrid_cmc:
         case integrator_t::lsf_mc:
            if(noise_type_set && coloured){
               err::warning("Warning - quantum:noise-type is ignored by Monte Carlo integrators.");
            }
            break;

         // the remaining spin dynamics integrators draw their own white noise;
         // the rotation integrators split it into two half steps per time
         // step, which a coloured process does not support
         default:
            if(coloured){
               err::fatal("Error - coloured quantum noise (quantum:noise-type = " + std::string(spectrum_name(noise_type)) +
                          ") is supported with sim:integrator = llg-heun, llg-quantum and spin-lattice only.");
            }
            break;

      }

      // the Ohmic bath belongs to the open-system thermostat
      if(mode != thermostat){
         for(size_t m = 0; m < internal::mp.size(); m++){
            if(internal::mp[m].alpha_e.get() > 0.0){
               err::warning("Warning - material:quantum-ohmic-damping is used by sim:integrator = llg-quantum only and is ignored here.");
               break;
            }
         }
      }

      if(mode == inactive) return;

      // GPU: the CUDA module runs the on-the-fly (HO) baths and the open-system
      // LLG on the device. Pre-generated noise, noise export, spin-lattice and
      // MPI builds stay on the CPU.
      if(gpu::acceleration){
         #ifdef MPICF
            err::fatal("Error - the quantum noise module on the GPU is not available in MPI builds; use the CUDA (serial) build.");
         #endif
         if(noise_generator == pre_generated){
            err::fatal("Error - quantum:noise-generator = pre-generated is not supported with gpu:acceleration; the GPU generates the noise on the fly (HO method).");
         }
         if(export_noise){
            err::fatal("Error - quantum:export-noise is not supported with gpu:acceleration; export the noise from a CPU run.");
         }
         if(mode == bath && integrate::integrator != integrate::integrator_t::llg_heun){
            err::fatal("Error - coloured quantum noise on the GPU is supported with sim:integrator = llg-heun and llg-quantum only.");
         }
      }

      // programs that generate their own thermal field bypass the bath
      if(mode == bath && (program::program == program::program_t::hamr_simulation ||
                          program::program == program::program_t::localised_temperature_pulse)){
         err::fatal("Error - coloured quantum noise is not supported with the hamr-simulation and "
                    "localised-temperature-pulse programs, which generate their own thermal field.");
      }

      // the open-system LLG integrates its own noise: switch the white thermal
      // field off so it draws no random numbers. A coloured bath keeps the
      // thermal field on; sim::calculate_thermal_fields() then copies the bath
      // sample into the thermal field arrays.
      if(mode == thermostat) sim::thermal_field_enabled = false;

      // a coloured bath uses the global temperature only
      if(coloured){
         bool rescaled = sim::local_temperature;
         for(int m = 0; m < mp::num_materials; m++){
            if(mp::material[m].temperature_rescaling_Tc > 0.0 && mp::material[m].temperature_rescaling_alpha != 1.0) rescaled = true;
         }
         if(rescaled){
            err::warning("Warning - temperature rescaling and material temperatures are not applied to coloured quantum noise; "
                         "sim:temperature is used.");
         }
      }

      if(export_noise) export_open((mode == thermostat) ? "auxiliary oscillator q" : "spin bath sample");

      std::ostringstream banner;
      banner << "Initialising quantum noise module: "
             << ((mode == thermostat) ? (ohmic ? "open-system LLG + Ohmic bath" : "open-system LLG") : "coloured bath") << ", spectrum "
             << spectrum_name(noise_type) << ", generator " << generator_name(noise_generator)
             << ", " << spin_bath.n_sites << " atoms";
      if(coloured && noise_generator == on_the_fly) banner << ", " << n_bath_modes << " bath modes";
      if(lattice_bath.n_sites > 0) banner << ", lattice bath";
      zlog << zTs() << banner.str() << std::endl;

      return;

   }

   //---------------------------------------------------------------------------
   // Function to release the FFTW resources of both baths
   //---------------------------------------------------------------------------
   void cleanup(){
      internal::fft_release(internal::spin_bath);
      internal::fft_release(internal::lattice_bath);
      return;
   }

} // end of quantum namespace
