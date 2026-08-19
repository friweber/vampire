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
//   Definitions for all internal module variables declared in internal.hpp.
//
//------------------------------------------------------------------------------

// Vampire headers
#include "quantum.hpp"

// Module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      //=====================================================================
      // Module state
      //=====================================================================

      bool enabled = false;

      std::vector<internal::mp_t> mp;

      noise_t noise_type = internal::quantum_zero;
      integration_t llg_method = internal::llg_fft;

      // Noise kind for llg-heun-quantum. Set via quantum:heun-noise-type.
      quantum::sld_noise::kind_t heun_noise_kind = quantum::sld_noise::kind_t::quantum;

      // Noise export parameters.
      // Default name is empty: a self-explanatory name like
      // "<noise-type>_<method>_noise.dat" is constructed in initialize()
      // if the user requested export without supplying an explicit filename.
      bool export_noise = false;
      std::string export_noise_filename = "";
      bool noise_export_header_written = false;
      uint64_t noise_export_burnin_steps = 0;
      uint64_t noise_export_step_count   = 0;

      //=====================================================================
      // Windowed noise generation state
      //=====================================================================

      bool windowed_mode = false;
      uint64_t window_size = 0;
      uint64_t M_decimation = 1;
      int window_start_fine = 0;
      int window_n_coarse = 0;
      int num_realizations = 0;
      double stored_dt_fine = 0.0;
      int stored_M = 1;
      double stored_T = 0.0;
      double noise_norm_factor = 0.0;
      double noise_inv_sqrt_S0 = 1.0;
      double noise_scale = 1.0;
      double noise_sigma = 0.0;
      int noise_seg = 0;
      int valid_n_coarse = 0;
      std::vector<double> prev_white_noise;
      std::vector<double> sqrt_PSD_window;

      //=====================================================================
      // Persistent FFT resources (used in windowed mode)
      //=====================================================================

      #ifdef FFT
      fftw_plan fft_forward = nullptr;
      fftw_plan fft_backward = nullptr;
      double* fft_in = nullptr;
      fftw_complex* fft_out = nullptr;
      double* fft_result = nullptr;
      #endif

      //=====================================================================
      // Per-material parameter arrays
      //=====================================================================

      std::vector<double> material_A_array;
      std::vector<double> material_gamma_array;
      std::vector<double> material_omega0_array;
      std::vector<double> material_S0_array;
      std::vector<double> material_inv_sqrt_S0_array;

      //=====================================================================
      // Noise field storage and indexing
      //=====================================================================

      std::vector<double> coarse_noise_field;

      std::vector<std::size_t> atom_idx_x;
      std::vector<std::size_t> atom_idx_y;
      std::vector<std::size_t> atom_idx_z;

      int noise_index = 0;

      //=====================================================================
      // Per-atom integration state arrays
      //=====================================================================

      std::vector<double> q_x_array;
      std::vector<double> q_y_array;
      std::vector<double> q_z_array;
      std::vector<double> p_x_array;
      std::vector<double> p_y_array;
      std::vector<double> p_z_array;

      std::vector<std::vector<double>> k1_storage;
      std::vector<std::vector<double>> k2_storage;
      std::vector<std::vector<double>> k3_storage;
      std::vector<std::vector<double>> k4_storage;
      std::vector<std::vector<double>> y_pred_storage;
      std::vector<std::vector<double>> y_in_storage;

      //=====================================================================
      // Auxiliary bath state (HO method, quantum / quantum-no-zero)
      // Single user-facing knob "quantum:bath-modes" — used for both the
      // Matsubara modes (orn-uhl) and the log-bath modes (log-bath-opt).
      //=====================================================================

      int n_bath_modes = 30;

      // Log-bath λ-range search tuning knobs (quantum_no_zero only).
      // Defaults match the cmp_noise reference; users can tighten / widen
      // via quantum:bath-scan-resolution / -omega-points / -decades.
      int    bath_scan_resolution   = 50;
      int    bath_scan_omega_points = 500;
      double bath_scan_decades      = 3.5;

      std::vector<double> matsubara_s_x;
      std::vector<double> matsubara_s_y;
      std::vector<double> matsubara_s_z;

      std::vector<double> qn_x_array;
      std::vector<double> qn_y_array;
      std::vector<double> qn_z_array;

      //=====================================================================
      // OU coefficients (quantum_zero noise type)
      // The PRNG is Vampire's mtrandom::grnd — declared in random.hpp,
      // seeded centrally by the parallel_rng_seed module.
      //=====================================================================

      double mats_white_amp = 0.0;

      std::vector<double> ou_decay;
      std::vector<double> ou_diffuse;
      std::vector<double> ou_drift;
      std::vector<double> ou_noise_amp;

      //=====================================================================
      // Opt-range log-bath state and coefficients (quantum_no_zero noise type)
      // (mode count lives above as n_bath_modes — same knob for both paths)
      //=====================================================================

      std::vector<double> lb_opt_s_x;
      std::vector<double> lb_opt_s_y;
      std::vector<double> lb_opt_s_z;
      std::vector<double> lb_opt_lambdas;
      std::vector<double> lb_opt_decay;
      std::vector<double> lb_opt_coeff_s;
      std::vector<double> lb_opt_amp;

      // Cache of the T_scaled value used to populate the coefficient arrays.
      // Set by setup_orn_uhl() / setup_log_bath_opt() immediately after they
      // fill the coefficients, so refresh_quantum_noise_for_T() can do a
      // straight equality compare without any sentinel logic.
      double last_T_scaled = 0.0;

      //=====================================================================
      // SLD phonon quantum noise — parallel OU bath for the lattice subsystem
      //=====================================================================

      std::vector<double> qn_phonon_x_array;
      std::vector<double> qn_phonon_y_array;
      std::vector<double> qn_phonon_z_array;

      std::vector<double> matsubara_phonon_s_x;
      std::vector<double> matsubara_phonon_s_y;
      std::vector<double> matsubara_phonon_s_z;

      std::vector<double> lb_opt_phonon_s_x;
      std::vector<double> lb_opt_phonon_s_y;
      std::vector<double> lb_opt_phonon_s_z;

      std::vector<double> ou_phonon_decay;
      std::vector<double> ou_phonon_diffuse;
      std::vector<double> ou_phonon_drift;
      std::vector<double> ou_phonon_noise_amp;
      double mats_phonon_white_amp = 0.0;

      std::vector<double> lb_opt_phonon_decay;
      std::vector<double> lb_opt_phonon_coeff_s;
      std::vector<double> lb_opt_phonon_amp;
      std::vector<double> lb_opt_phonon_lambdas;

      std::vector<double> material_phonon_amp_array;
      std::vector<double> material_phonon_omega0_array;
      std::vector<double> material_phonon_gamma_array;

      std::vector<double> material_sld_spin_amp_array;
      std::vector<double> material_sld_classical_sigma_array;

      int n_cascade_modes = 1;
      std::vector<double> lb_opt_casc_s_x;
      std::vector<double> lb_opt_casc_s_y;
      std::vector<double> lb_opt_casc_s_z;
      std::vector<double> lb_opt_casc_coeff;

      // Butterworth LP post-filter state (quantum_no_zero filtered path)
      BiquadCoeffs butter_stage[2]    = {};
      bool         butter_coeffs_valid = false;
      double       butter_cutoff_Hz    = 0.0;
      std::vector<double> filter_state_x;
      std::vector<double> filter_state_y;
      std::vector<double> filter_state_z;

      double last_T_scaled_phonon = 0.0;

      //=====================================================================
      // SLD / ASD FFT spin noise (pre-generated, pure coth / coth-1 shape)
      //=====================================================================

      std::vector<double> sld_fft_spin_x;
      std::vector<double> sld_fft_spin_y;
      std::vector<double> sld_fft_spin_z;
      uint64_t sld_fft_step_index = 0;
      uint64_t sld_fft_n_fine     = 0;
      int      sld_fft_num_atoms  = 0;

      //=====================================================================
      // SLD / ASD spin-noise export state
      //=====================================================================

      bool        sld_export_noise          = false;
      std::string sld_export_noise_filename = "";
      int         sld_export_noise_atom     = 0;

   } // end of internal namespace

} // end of quantum namespace
