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

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
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
         return constants::kB/(constants::hbar*constants::gyromagnetic_ratio);
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

         switch(program::program){
            case  0: n_steps = tt;      return true; // benchmark
            case  1: n_steps = et + tt; return true; // time series
            case  2: n_steps = et + lt; return true; // hysteresis
            case  3: n_steps = et + lt; return true; // static hysteresis
            case  4: n_steps = et + lt; return true; // Curie temperature
            case  5: n_steps = et + tt; return true; // field cool
            case  6: n_steps = et + tt; return true; // laser / temperature pulse
            case  7: n_steps = et + tt; return true; // HAMR
            case 11: n_steps = et + tt; return true; // Lagrange multiplier
            case 12: n_steps = et + lt; return true; // partial hysteresis
            case 13: n_steps = et + tt; return true; // localised temperature pulse
            case 14: n_steps = et + tt; return true; // effective damping
            case 15: n_steps = et + tt; return true; // FMR
            case 16: n_steps = et + tt; return true; // local field cool
            case 17: n_steps = et + tt; return true; // electrical pulse
            case 18: n_steps = et + tt; return true; // field pulse
            case 52: n_steps = et + tt; return true; // domain walls
            case 70: n_steps = et + lt; return true; // field sweep
            case 74: n_steps = et + tt; return true; // spin waves
            default: return false;
         }

      }

      //---------------------------------------------------------------------------
      // Function to report whether the program moves sim::temperature during
      // the run; pre-generated noise is shaped once and cannot follow it
      //---------------------------------------------------------------------------
      bool dynamic_temperature_program(){

         switch(program::program){
            case  5: // field cool
            case  6: // laser / temperature pulse
            case  7: // HAMR
            case 13: // localised temperature pulse
            case 16: // local field cool
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
            terminaltextcolor(RED);
            std::cerr << "Error: program " << program::program
                      << " is not supported with pre-generated quantum noise." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         if(dynamic_temperature_program()){
            terminaltextcolor(RED);
            std::cerr << "Error: pre-generated quantum noise is shaped once, at the starting temperature, "
                      << "and cannot follow a program that changes the temperature. "
                      << "Use quantum:noise-generator = on-the-fly." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         return n_run;

      }

      //---------------------------------------------------------------------------
      // Thermostat (llg-quantum): per-material Lorentzian response and a spin
      // bath whose amplitude balances that response. Classical noise is
      // allowed here; it is the classical open-system model.
      //---------------------------------------------------------------------------
      static void initialize_thermostat(){

         const int n_mat = mp::num_materials;
         if(mp.size() < static_cast<size_t>(n_mat)) mp.resize(n_mat);

         material_A.resize(n_mat);
         material_gamma.resize(n_mat);
         material_omega0.resize(n_mat);
         material_S0.resize(n_mat);

         std::vector<double> amp(n_mat);

         for(int m = 0; m < n_mat; m++){

            const double alpha = mp::material[m].alpha;
            const double gamma = mp[m].gamma.get();
            const double omega0 = mp[m].omega0.get();

            if(gamma <= 0.0 || omega0 <= 0.0){
               terminaltextcolor(RED);
               std::cerr << "Error: material " << m + 1 << " needs quantum-lorentzian-width and "
                         << "quantum-lorentzian-central-frequency > 0 for sim:integrator = llg-quantum." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }

            // the moment in units of hbar*gamma_e converts the reduced temperature to a field
            const double S0 = mp::material[m].mu_s_SI/(constants::hbar*constants::gyromagnetic_ratio);
            const double A = alpha*std::pow(omega0, 4)/gamma;

            material_A[m] = A;
            material_gamma[m] = gamma;
            material_omega0[m] = omega0;
            material_S0[m] = S0;

            // on-the-fly: the oscillator filters the bath, amplitude sqrt(gamma A/S0);
            // pre-generated: the filter is in the spectrum, amplitude 1/sqrt(S0)
            amp[m] = (noise_generator == pre_generated) ? 1.0/std::sqrt(S0) : std::sqrt(gamma*A/S0);

         }

         const uint64_t n_run = pre_generated_run_length();
         const int n = num_local_atoms();

         setup_bath(spin_bath, noise_type, noise_generator, n, mp::dt, spin_temperature_scale(), amp, amp, true, n_run);
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
         const double pref = std::sqrt(mp::dt*constants::hbar*constants::gyromagnetic_ratio/(2.0*constants::kB));

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
            terminaltextcolor(RED);
            std::cerr << "Error: quantum noise for spin-lattice dynamics requires the spin-lattice module "
                      << "to be initialised first." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         const int n_mat = mp::num_materials;
         const int n = num_local_atoms();

         // spin bath
         const double pref = std::sqrt(mp::dt*constants::hbar*constants::gyromagnetic_ratio/(2.0*constants::kB));
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

      mode = inactive;

      // sim:integrator = llg-heun-quantum is the old spelling of llg-heun with a
      // quantum bath; it implied the zero-point spectrum
      if(sim::integrator == sim::llg_heun_quantum){
         zlog << zTs() << "Warning: sim:integrator = llg-heun-quantum is deprecated; use llg-heun with quantum:noise-type." << std::endl;
         if(!noise_type_set) noise_type = quantum_zero;
      }

      switch(sim::integrator){

         case sim::llg_quantum:
            initialize_thermostat();
            break;

         case sim::llg_heun:
         case sim::llg_heun_quantum:
            if(noise_type != classical) initialize_llg_bath();
            break;

         case sim::suzuki_trotter:
            if(noise_type != classical) initialize_lattice_baths();
            break;

         default:
            if(noise_type_set && noise_type != classical){
               zlog << zTs() << "Warning: quantum:noise-type is ignored by the selected integrator." << std::endl;
            }
            break;

      }

      if(mode == inactive) return;

      // the module provides the thermal noise from here on
      sim::hamiltonian_simulation_flags[3] = 0;

      // a coloured bath uses the global temperature only
      if(noise_type != classical){
         bool rescaled = sim::local_temperature;
         for(int m = 0; m < mp::num_materials; m++){
            if(mp::material[m].temperature_rescaling_Tc > 0.0 && mp::material[m].temperature_rescaling_alpha != 1.0) rescaled = true;
         }
         if(rescaled){
            zlog << zTs() << "Warning: temperature rescaling and material temperatures are not applied to coloured quantum noise; "
                 << "sim:temperature is used." << std::endl;
         }
      }

      if(export_noise) export_open((mode == thermostat) ? "auxiliary oscillator q" : "spin bath sample");

      zlog << zTs() << "Initialising quantum noise module: "
           << ((mode == thermostat) ? "open-system LLG" : "coloured bath") << ", spectrum "
           << spectrum_name(noise_type) << ", generator " << generator_name(noise_generator)
           << ", " << spin_bath.n_sites << " atoms";
      if(noise_type != classical && noise_generator == on_the_fly) zlog << ", " << n_bath_modes << " bath modes";
      if(lattice_bath.n_sites > 0) zlog << ", lattice bath";
      zlog << std::endl;

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
