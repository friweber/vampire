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
//   Initialization for the quantum thermostat module.
//
//   initialize()     — Common entry point: populates material parameter arrays,
//                      disables the standard thermal field, and dispatches to
//                      the method-specific initializer (FFT or HO).
//
//   initialize_FFT() — Sets up the FFT noise pipeline. Determines whether to
//                      use windowed or non-windowed mode based on window_size
//                      vs total simulation length. Allocates all integration
//                      arrays and generates the initial noise.
//
//   initialize_HO()  — Allocates integration arrays for the on-the-fly
//                      harmonic oscillator noise method.
//
//   supported_program() — Returns true if the selected program is compatible
//                         with the quantum thermostat, and outputs the total
//                         number of simulation time steps.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "quantum.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// Module headers
#include "internal.hpp"

namespace quantum{

   // Forward declarations for method-specific initialization
   bool supported_program(uint64_t& total_simulation_time);
   void initialize_FFT();
   void initialize_HO();

   //========================================================================
   // Common initialization entry point
   //========================================================================
   void initialize(){

      using namespace internal;

      if(!enabled) return;

      // Disable standard thermal field — quantum module provides its own noise
      sim::hamiltonian_simulation_flags[3] = 0;

      //---------------------------------------------------------------------
      // Populate per-material parameter arrays from input file values
      //---------------------------------------------------------------------
      for(int m = 0; m < mp::num_materials; m++){

         double alpha  = mp::material[m].alpha;
         double gamma  = internal::mp[m].gamma.get();
         double omega0 = internal::mp[m].omega0.get();

         // Validate: gamma and omega0 must be strictly positive. Silent zeros
         // produce A = NaN/0 and divergent noise downstream (FDT prefactor
         // blows up, coth(ω/0) is ill-defined).
         if (gamma <= 0.0 || omega0 <= 0.0) {
            std::cerr << "Error: quantum thermostat material " << m
                      << " has invalid Lorentzian parameters:\n"
                      << "  quantum-lorentzian-width            (Gamma)  = " << gamma  << "\n"
                      << "  quantum-lorentzian-central-frequency (omega0) = " << omega0 << "\n"
                      << "Both must be strictly positive (rad/s). Check the .mat file."
                      << std::endl;
            err::vexit();
         }

         // Spin magnitude: S0 = mu_s / mu_B
         double S0 = mp::material[m].mu_s_SI / 9.274009994e-24;
         double inv_sqrt_S0 = (S0 > 0.0) ? 1.0 / std::sqrt(S0) : 1.0;

         // Lorentzian amplitude: A = alpha * omega0^4 / Gamma
         double A = alpha * pow(omega0, 4) / gamma;

         material_gamma_array.push_back(gamma);
         material_omega0_array.push_back(omega0);
         material_A_array.push_back(A);
         material_S0_array.push_back(S0);
         material_inv_sqrt_S0_array.push_back(inv_sqrt_S0);
      }

      //---------------------------------------------------------------------
      // Construct default export-noise filename if the user requested
      // export but did not supply a name. Self-explanatory format:
      //   <noise-type>_<method>_noise.dat
      //---------------------------------------------------------------------
      if (export_noise && export_noise_filename.empty()) {
         export_noise_filename = std::string(noise_type_name(noise_type)) + "_"
                               + llg_method_short_name(llg_method) + "_noise.dat";
      }

      //---------------------------------------------------------------------
      // Unified startup banner. Same shape for FFT and HO; method-specific
      // lines are printed by the per-method initializer.
      //---------------------------------------------------------------------
      std::cout << "Initializing Quantum Noise Module..." << std::endl;
      std::cout << "  Method            : " << llg_method_long_name(llg_method) << std::endl;
      std::cout << "  Noise type        : " << noise_type_name(noise_type)      << std::endl;
      std::cout << "  Temperature       : " << sim::temperature                 << " K"  << std::endl;
      std::cout << "  Time step         : " << mp::dt                           << " s"  << std::endl;
      if (export_noise) {
         std::cout << "  Noise export      : " << export_noise_filename << std::endl;
      }
      std::cout << "  Materials         : " << mp::num_materials << std::endl;
      for (int m = 0; m < mp::num_materials; m++) {
         std::cout << "    [" << m << "]"
                   << " A="  << material_A_array[m]
                   << "  Gamma=" << material_gamma_array[m]
                   << "  omega0=" << material_omega0_array[m]
                   << "  S0=" << material_S0_array[m] << std::endl;
      }

      //---------------------------------------------------------------------
      // Dispatch to method-specific initialization
      //---------------------------------------------------------------------
      if(llg_method == llg_ho){
         initialize_HO();
      }
      else if(llg_method == llg_fft){
         initialize_FFT();
      }

      return;
   }

   //========================================================================
   // FFT-specific initialization
   //
   // Determines the total simulation length, computes coarse grid parameters,
   // allocates integration arrays, and generates noise. If the requested
   // window_size is smaller than the total simulation, windowed mode is used.
   //========================================================================
   void initialize_FFT(){

      using namespace internal;

      // Check that this program is supported by the quantum thermostat
      uint64_t total_simulation_time = 0;
      if(!supported_program(total_simulation_time)){
         std::cerr << "Error: program " << program::program
                   << " is not supported with the quantum thermostat." << std::endl;
         err::vexit();
      }

      // Convert temperature from Kelvin to internal units (rad/s).
      const double T = scale_temperature(sim::temperature);

      // If window_size is zero (default), set to full simulation length
      if(window_size == 0){
         window_size = total_simulation_time + 1;
      }

      // Fine time grid (used for spin dynamics integration)
      double dt_fine = mp::dt;
      uint64_t n_fine = total_simulation_time + 1;

      // Coarse time grid (used for noise generation)
      int M = internal::M_decimation;
      int n_coarse = (n_fine > 0) ? ((n_fine - 1) / M + 1) : 0;

      // FFT-specific banner block (matches the unified banner shape)
      std::cout << "  Total time steps  : " << n_fine    << std::endl;
      std::cout << "  Window size       : " << window_size << "  (fine steps)" << std::endl;
      std::cout << "  Interpolation M   : " << M         << "  (-> " << n_coarse << " coarse steps)" << std::endl;
      std::cout << "  Note: FFT noise is fixed at the initial sim::temperature." << std::endl;
      std::cout << "        Use quantum:llg-method=llg-ho for dynamic-T programs." << std::endl;

      // Number of atoms (including MPI boundary atoms)
      #ifdef MPICF
         const int num_atoms_total = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      #else
         const int num_atoms_total = atoms::num_atoms;
      #endif

      int realizations = num_atoms_total * 3;  // 3 spatial components per atom

      //---------------------------------------------------------------------
      // Allocate per-atom integration arrays (9 components: S, q, p)
      //---------------------------------------------------------------------
      q_x_array.resize(num_atoms_total, 0.0);
      q_y_array.resize(num_atoms_total, 0.0);
      q_z_array.resize(num_atoms_total, 0.0);
      p_x_array.resize(num_atoms_total, 0.0);
      p_y_array.resize(num_atoms_total, 0.0);
      p_z_array.resize(num_atoms_total, 0.0);

      k1_storage.resize(num_atoms_total, std::vector<double>(9));
      k2_storage.resize(num_atoms_total, std::vector<double>(9));
      k3_storage.resize(num_atoms_total, std::vector<double>(9));
      k4_storage.resize(num_atoms_total, std::vector<double>(9));
      y_pred_storage.resize(num_atoms_total, std::vector<double>(9));
      y_in_storage.resize(num_atoms_total, std::vector<double>(9));

      // Per-atom noise carrier — written once per RK4 step by
      // draw_noise_all_atoms_FFT(), reused across K1-K4 via collect_H_FFT().
      // (Same array the HO path uses — single noise carrier across both methods.)
      qn_x_array.resize(num_atoms_total, 0.0);
      qn_y_array.resize(num_atoms_total, 0.0);
      qn_z_array.resize(num_atoms_total, 0.0);

      //---------------------------------------------------------------------
      // Choose between windowed and non-windowed noise generation
      //---------------------------------------------------------------------
      if(window_size < n_fine){

         // --- Windowed mode with overlap-save ---
         windowed_mode = true;
         int window_n_coarse_local = (window_size > 0) ? ((window_size - 1) / M + 1) : 0;

         // Round up to the next multiple of OVERLAP_SAVE_SEGMENTS so the
         // window divides cleanly into segments for the overlap-save scheme.
         if (window_n_coarse_local % OVERLAP_SAVE_SEGMENTS != 0) {
            window_n_coarse_local = ((window_n_coarse_local / OVERLAP_SAVE_SEGMENTS) + 1)
                                    * OVERLAP_SAVE_SEGMENTS;
         }
         const int seg      = window_n_coarse_local / OVERLAP_SAVE_SEGMENTS;
         const int valid_nc = OVERLAP_SAVE_VALID * seg;

         std::cout << "  Windowed noise generation enabled (overlap-save):" << std::endl;
         std::cout << "    Window fine steps: " << window_size << std::endl;
         std::cout << "    FFT coarse steps: " << window_n_coarse_local
                   << " (" << OVERLAP_SAVE_SEGMENTS << " segments of " << seg << ")" << std::endl;
         std::cout << "    Valid coarse steps per window: " << valid_nc << std::endl;

         assign_unique_indices(valid_nc, num_atoms_total);
         init_noise_structures(window_n_coarse_local, realizations, dt_fine, M, T);
         generate_noise_window();
      }
      else{

         // --- Non-windowed mode (all noise generated at once) ---
         windowed_mode = false;
         window_start_fine = 0;

         assign_unique_indices(n_coarse, num_atoms_total);
         init_noise_structures(n_coarse, realizations, dt_fine, M, T);
         if (n_fine > 0) {
            calculate_noise(realizations, dt_fine, M, T, n_coarse, coarse_noise_field);
         }
      }

      // Export noise if requested (non-windowed mode only;
      // windowed mode exports automatically after each window generation)
      if(export_noise && !windowed_mode){
         #ifdef MPICF
         if(vmpi::my_rank == 0){
         #endif
            std::cout << "  Exporting noise data for analysis..." << std::endl;
            export_noise_data(export_noise_filename);
         #ifdef MPICF
         }
         #endif
      }

      return;
   }

   //========================================================================
   // Harmonic Oscillator-specific initialization
   //
   // Only allocates integration arrays. HO noise is generated on-the-fly
   // during integration (no pre-generation needed).
   //========================================================================
   void initialize_HO(){

      using namespace internal;

      // HO-specific banner block (matches the unified banner shape).
      // HO does not pre-generate a noise field, so the "total time steps"
      // line is just the user's requested run length for reference.
      const uint64_t total_steps = sim::equilibration_time + sim::total_time;
      std::cout << "  Total time steps  : " << total_steps << std::endl;
      if (noise_type != classical) {
         std::cout << "  Bath modes        : " << n_bath_modes << std::endl;
      }

      // HO noise-export burn-in.  The auxiliary oscillator and bath modes
      // start at zero; the first ~few hundred steps are a pure transient
      // that, if recorded, dominates the low-frequency end of the Welch
      // PSD and destroys the shape match.  Drop min(total_steps/5, 20000)
      // samples before recording — matches cmp_noise's convention.
      if (export_noise) {
         noise_export_burnin_steps = std::min<uint64_t>(total_steps / 5, 20000);
         noise_export_step_count   = 0;
         std::cout << "  Export burn-in    : " << noise_export_burnin_steps
                   << " steps (HO transient skipped)" << std::endl;
      }

      #ifdef MPICF
         const int num_atoms_total = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      #else
         const int num_atoms_total = atoms::num_atoms;
      #endif

      q_x_array.resize(num_atoms_total, 0.0);
      q_y_array.resize(num_atoms_total, 0.0);
      q_z_array.resize(num_atoms_total, 0.0);
      p_x_array.resize(num_atoms_total, 0.0);
      p_y_array.resize(num_atoms_total, 0.0);
      p_z_array.resize(num_atoms_total, 0.0);

      k1_storage.resize(num_atoms_total, std::vector<double>(9));
      k2_storage.resize(num_atoms_total, std::vector<double>(9));
      k3_storage.resize(num_atoms_total, std::vector<double>(9));
      k4_storage.resize(num_atoms_total, std::vector<double>(9));
      y_pred_storage.resize(num_atoms_total, std::vector<double>(9));
      y_in_storage.resize(num_atoms_total, std::vector<double>(9));

      // Allocate per-atom noise arrays (reused across K1-K4 each step)
      qn_x_array.resize(num_atoms_total, 0.0);
      qn_y_array.resize(num_atoms_total, 0.0);
      qn_z_array.resize(num_atoms_total, 0.0);

      // Dispatch to the noise-type-specific setup (allocates per-atom
      // oscillator state and precomputes coefficients). Lives in
      // noise_ho.cpp alongside the matching generator.
      // (mtrandom::grnd is already seeded by Vampire's parallel_rng_seed
      //  module — no extra seeding is needed here.)
      if (noise_type != classical) {
         if (noise_type == quantum_no_zero) {
            // When a Butterworth LP cutoff is set, use the filtered setup
            // (which calls setup_log_bath_opt internally then allocates filter state).
            if (butter_cutoff_Hz > 0.0)
               setup_log_bath_filtered(num_atoms_total, butter_cutoff_Hz / 1.0e12);
            else
               setup_log_bath_opt(num_atoms_total);
         } else {
            setup_orn_uhl(num_atoms_total);
         }
      }

      return;
   }

   //========================================================================
   // Programs that move sim::temperature while the simulation runs.
   //========================================================================
   bool internal::dynamic_temperature_program(){
      switch (program::program) {
         case  5:   // field cool
         case  6:   // laser / temperature pulse
         case  7:   // HAMR
         case 13:   // localised temperature pulse
         case 16:   // local field cool
            return true;
         default:
            return false;
      }
   }

   //========================================================================
   // Check if the selected simulation program is supported by the quantum
   // thermostat module's currently selected llg_method, and report the
   // total simulation time-step count.
   //
   // Two classes of programs:
   //   * Static-T or single-pulse: supported on both FFT and HO paths.
   //   * Dynamic-T (field cool, HAMR, localised T pulse, local field cool):
   //     supported only on HO, because refresh_quantum_noise_for_T() tracks
   //     sim::temperature every step. The FFT noise field is pre-generated
   //     at init T and cannot follow runtime T changes, so it stays gated.
   //========================================================================
   bool supported_program(uint64_t& total_simulation_time){

      const uint64_t et = sim::equilibration_time;
      const uint64_t tt = sim::total_time;
      const uint64_t lt = sim::loop_time;
      const bool ho     = (internal::llg_method == internal::llg_ho);

      switch (program::program) {

         // --- Static-T or single-pulse programs (both paths) ---
         case  0: total_simulation_time = tt;      return true;  // benchmark
         case  1: total_simulation_time = et + tt; return true;  // time series
         case  2: total_simulation_time = et + lt; return true;  // hysteresis
         case  3: total_simulation_time = et + lt; return true;  // static hysteresis
         case  4: total_simulation_time = et + lt; return true;  // Curie temperature
         case 11: total_simulation_time = et + tt; return true;  // LaGrange multiplier
         case 12: total_simulation_time = et + lt; return true;  // partial hysteresis
         case 14: total_simulation_time = et + tt; return true;  // effective damping
         case 15: total_simulation_time = et + tt; return true;  // FMR
         case 17: total_simulation_time = et + tt; return true;  // electrical pulse
         case 18: total_simulation_time = et + tt; return true;  // field pulse
         case 52: total_simulation_time = et + tt; return true;  // domain walls
         case 70: total_simulation_time = et + lt; return true;  // field sweep
         case 74: total_simulation_time = et + tt; return true;  // spin waves

         // --- Dynamic-T programs: HO only (FFT pipeline can't track runtime T) ---
         case  5: total_simulation_time = et + tt; return ho;    // field cool
         case  6: total_simulation_time = et + tt; return ho;    // laser / temperature pulse
         case  7: total_simulation_time = et + tt; return ho;    // HAMR
         case 13: total_simulation_time = et + tt; return ho;    // localised T pulse
         case 16: total_simulation_time = et + tt; return ho;    // local field cool

         default:                                  return false;
      }
   }

} // end of quantum namespace
