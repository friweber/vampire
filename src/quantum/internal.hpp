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
//   Quantum Thermostat Module — Internal Header
//
//   This module implements a quantum thermostat for atomistic spin dynamics
//   based on colored noise with a Lorentzian spectral density. The thermal
//   noise field is generated via one of two methods:
//
//     FFT method:  Pre-generate noise on a coarse time grid using FFT-based
//                  spectral shaping, then interpolate to the fine integration
//                  time step. Supports windowed generation for long runs.
//
//     HO method:   On-the-fly colored noise via auxiliary oscillator modes.
//                  Two noise types are supported:
//                    quantum_zero    -> exact Ornstein-Uhlenbeck (n_bath_modes Matsubara modes)
//                    quantum_no_zero -> log-spaced OU bath (n_bath_modes modes,
//                                      lambda range found at init by 2-D grid scan)
//
//   Both methods use a 9-component RK4 integrator per atom:
//     y = (S_x, S_y, S_z, q_x, q_y, q_z, p_x, p_y, p_z)
//   where S is the spin unit vector, q is an auxiliary position (oscillator
//   displacement), and p is an auxiliary momentum.
//
//   The power spectral density (PSD) of the noise is:
//
//     Classical:       P(w) = 2T * A * G / [(w0^2 - w^2)^2 + G^2 * w^2]
//     Quantum:         P(w) = coth(w/2T) * A * G * w / [(w0^2 - w^2)^2 + G^2 * w^2]
//     Quantum-no-zero: P(w) = [coth(w/2T) - 1] * A * G * w / [(w0^2 - w^2)^2 + G^2 * w^2]
//
//   where A = alpha * w0^4 / G, with alpha the Gilbert damping constant.
//
//   References:
//     Barker & Bauer, Phys. Rev. B 100, 140401(R) (2019)
//     https://doi.org/10.21105/joss.06263
//
//------------------------------------------------------------------------------

#ifndef QUANTUM_INTERNAL_H_
#define QUANTUM_INTERNAL_H_

// C++ standard library headers
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#ifdef FFT
#include <fftw3.h>
#endif

// Vampire headers
#include "constants.hpp"
#include "quantum.hpp"

namespace quantum{

   namespace internal{

      //=====================================================================
      // Temperature scaling — single, authoritative conversion.
      //
      // Converts a temperature in Kelvin to the internal "energy/ℏ" units
      // (rad/s) that ω₀, Γ, and all noise PSDs are expressed in:
      //
      //     T_scaled = T_K · kB / (ℏ · γ)
      //
      // ALL call sites that need a scaled temperature must call this helper
      // — never cache the result.  Programs that vary sim::temperature at
      // runtime (temperature pulse, laser pulse, …) get correct dynamic-T
      // behaviour only because of the per-call reading of sim::temperature.
      //=====================================================================
      inline double scale_temperature(const double T_K) {
         return T_K * constants::kB / (constants::hbar * constants::gyromagnetic_ratio);
      }

      //=====================================================================
      // Phonon temperature scaling — returns thermal frequency in rad/ps.
      //
      // Used by the SLD phonon quantum noise bath so that its Matsubara
      // frequencies are on the same time scale as the phonon integrator
      // (which uses dt_ps = mp::dt_SI * 1e12).
      //
      //     T_phonon = T_K * kB_eV / hbar_eVps
      //
      // with hbar_eVps = hbar [J*s] * 1e12 / e_charge [J/eV] ≈ 6.582e-4 eV*ps
      //=====================================================================
      inline double scale_temperature_phonon(const double T_K) {
         static const double hbar_eVps = constants::hbar * 1e12 / 1.6021766208e-19;
         return T_K * constants::kB_eV / hbar_eVps;
      }

      //=====================================================================
      // Internal data type definitions
      //=====================================================================

      //---------------------------------------------------------------------
      // Simple wrapper class for a double variable with a "set" flag.
      // Used for material parameters that may or may not be specified
      // in the input file.
      //---------------------------------------------------------------------
      class set_double_t{

      private:
         double value;
         bool setf;

      public:
         set_double_t() : value(0.0), setf(false) {}

         void set(double in_value){ value = in_value; setf = true; }
         double get(){ return value; }
         bool is_set(){ return setf; }
      };

      //---------------------------------------------------------------------
      // Material parameters class for quantum thermostat.
      // Stores the Lorentzian width (Gamma) and central frequency (omega0)
      // per material. These are read from the material file.
      //---------------------------------------------------------------------
      class mp_t{

      public:
         set_double_t gamma;   // Lorentzian width parameter Gamma [rad/s]
         set_double_t omega0;  // Lorentzian central frequency omega0 [rad/s]

         mp_t(const unsigned int max_materials = 100){
            gamma.set(0.0);
            omega0.set(0.0);
         }
      };

      //=====================================================================
      // Enumerations
      //=====================================================================

      /// Noise spectral density type
      enum noise_t {
         classical,        // Classical white noise (PSD proportional to T)
         quantum_zero,     // Quantum noise including zero-point fluctuations
         quantum_no_zero   // Quantum noise excluding zero-point term
      };

      /// Stringified noise type — for startup banner and default file names
      inline const char* noise_type_name(noise_t nt) {
         switch (nt) {
            case classical:       return "classical";
            case quantum_zero:    return "quantum";
            case quantum_no_zero: return "quantum-no-zero";
         }
         return "unknown";
      }

      /// LLG integration method
      enum integration_t {
         llg_fft,   // FFT-based pre-generated colored noise
         llg_ho     // Harmonic oscillator on-the-fly colored noise
      };

      /// Stringified method — lowercase for default file names, uppercase for banner
      inline const char* llg_method_short_name(integration_t m) {
         return m == llg_ho ? "ho" : "fft";
      }
      inline const char* llg_method_long_name(integration_t m) {
         return m == llg_ho ? "HO  (on-the-fly noise generation)"
                            : "FFT (pre-generated noise field)";
      }

      //=====================================================================
      // Module state variables
      //=====================================================================

      extern bool enabled;                        // Module activation flag
      extern std::vector<internal::mp_t> mp;      // Per-material input parameters

      extern noise_t noise_type;                  // Selected noise spectral type
      extern integration_t llg_method;            // Selected integration method

      // Noise kind for the llg-heun-quantum integrator.
      // Set via quantum:heun-noise-type — does NOT enable the full quantum thermostat.
      // Default: quantum (coth, HO on-the-fly).
      extern quantum::sld_noise::kind_t heun_noise_kind;
      extern uint64_t window_size;                // Noise window size in fine steps (0 = full run)
      extern uint64_t M_decimation;               // Coarse-to-fine interpolation factor

      extern bool export_noise;                   // Flag to export noise data to file
      extern std::string export_noise_filename;   // Filename for noise export
      extern bool noise_export_header_written;     // True after header has been written

      // HO-export burn-in: the auxiliary oscillator (q, p) and the bath modes s_k
      // all start at zero, so the first ~few hundred steps are a pure transient.
      // export_ho_noise_step() drops the first noise_export_burnin_steps calls
      // before recording, matching cmp_noise's min(nsteps/5, 20000) convention.
      // noise_export_step_count tracks how many times export_ho_noise_step has
      // been called (incremented every step regardless of whether anything was
      // written).  No effect on the FFT path, which generates a stationary
      // pre-baked field with no warm-up.
      extern uint64_t noise_export_burnin_steps;
      extern uint64_t noise_export_step_count;

      //=====================================================================
      // FFT resources (persistent in windowed mode)
      //=====================================================================

      #ifdef FFT
      extern fftw_plan fft_forward;     // Forward  FFT plan (real -> complex)
      extern fftw_plan fft_backward;    // Inverse  FFT plan (complex -> real)
      extern double* fft_in;            // Input  array for forward FFT
      extern fftw_complex* fft_out;     // Output array for forward FFT
      extern double* fft_result;        // Output array for inverse FFT
      #endif

      //=====================================================================
      // Overlap-save bookkeeping constants (windowed FFT noise generation)
      //
      // The FFT-domain window is divided into OVERLAP_SAVE_SEGMENTS equal
      // segments. After each forward+inverse FFT, the first and last
      // segments are corrupted by circular-convolution wrap-around and
      // are discarded; only the middle OVERLAP_SAVE_VALID segments are
      // kept as valid coarse-grid noise samples. The fixed 33 % overhead
      // (2 wasted segments out of 6) is the price of overlap-save and is
      // independent of window size or PSD shape.
      //
      //   window_n_coarse = OVERLAP_SAVE_SEGMENTS · noise_seg
      //   valid_n_coarse  = OVERLAP_SAVE_VALID    · noise_seg
      //
      // OVERLAP_SAVE_SEGMENTS must divide window_n_coarse evenly.  The
      // input parser enforces (window_size % OVERLAP_SAVE_SEGMENTS) == 0.
      //=====================================================================
      constexpr int OVERLAP_SAVE_SEGMENTS = 6;
      constexpr int OVERLAP_SAVE_VALID    = 4;

      //=====================================================================
      // Windowed noise generation state
      //
      // When the simulation is longer than the window size, noise is
      // generated in successive windows rather than all at once. This
      // keeps memory usage bounded for long simulations.
      //=====================================================================

      extern bool windowed_mode;                  // True if using windowed noise generation
      extern int window_start_fine;               // Fine step index where current window starts
      extern int window_n_coarse;                 // Number of coarse steps per window
      extern int num_realizations;                // Number of noise realizations (3 * num_atoms)

      // Cached parameters for window regeneration
      extern double stored_dt_fine;               // Fine time step dt [s]
      extern int stored_M;                        // Decimation factor M
      extern double stored_T;                     // Scaled temperature T [energy units]

      // Precomputed scaling constants (set once in init_noise_structures)
      extern double noise_norm_factor;            // 1/n_coarse — FFT normalization
      extern double noise_inv_sqrt_S0;            // 1/sqrt(S0) — spin magnitude normalization
      extern double noise_scale;                  // sqrt(dt_fine/dt_coarse) — variance matching
      extern double noise_sigma;                  // 1/sqrt(dt_fine) — white noise std dev
      extern int noise_seg;                       // Overlap-save segment size (window_n_coarse / 6)
      extern int valid_n_coarse;                  // Valid coarse samples per window (4 * noise_seg)
      extern std::vector<double> prev_white_noise; // White noise tail from previous window (2*seg per realization)
      extern std::vector<double> sqrt_PSD_window; // Precomputed sqrt(PSD) for window frequencies

      //=====================================================================
      // Per-material parameter arrays (populated during initialization)
      //=====================================================================

      extern std::vector<double> material_A_array;              // Lorentzian amplitude A
      extern std::vector<double> material_gamma_array;          // Lorentzian width Gamma [rad/s]
      extern std::vector<double> material_omega0_array;         // Lorentzian central frequency omega0 [rad/s]
      extern std::vector<double> material_S0_array;             // Spin magnitude S0 = mu_s / mu_B
      extern std::vector<double> material_inv_sqrt_S0_array;    // 1/sqrt(S0)

      //=====================================================================
      // Noise field storage
      //=====================================================================

      extern std::vector<double> coarse_noise_field;   // Flat 1D array [realization * n_coarse]

      // Per-atom index offsets into coarse_noise_field for each spatial component
      extern std::vector<std::size_t> atom_idx_x;
      extern std::vector<std::size_t> atom_idx_y;
      extern std::vector<std::size_t> atom_idx_z;

      extern int noise_index;   // Current fine time step index (advances each integration step)

      //=====================================================================
      // Per-atom integration state arrays
      //=====================================================================

      // Auxiliary oscillator position q and momentum p
      extern std::vector<double> q_x_array, q_y_array, q_z_array;
      extern std::vector<double> p_x_array, p_y_array, p_z_array;

      // RK4 temporary storage: k1-k4 slopes, predicted state, and initial state
      // Each inner vector has 9 components: (S_x, S_y, S_z, q_x, q_y, q_z, p_x, p_y, p_z)
      extern std::vector<std::vector<double>> k1_storage;
      extern std::vector<std::vector<double>> k2_storage;
      extern std::vector<std::vector<double>> k3_storage;
      extern std::vector<std::vector<double>> k4_storage;
      extern std::vector<std::vector<double>> y_pred_storage;
      extern std::vector<std::vector<double>> y_in_storage;

      //=====================================================================
      // Auxiliary bath state for the HO method (quantum / quantum-no-zero).
      //
      // n_bath_modes is the single user-facing knob for both:
      //   quantum         -> n_bath_modes Matsubara modes (orn-uhl)
      //   quantum-no-zero -> n_bath_modes log-spaced OU modes (log-bath-opt)
      //
      // Flat per-atom state arrays: size = num_atoms * n_bath_modes
      // Index: atom * n_bath_modes + mode
      // Used by generate_quantum_noise_orn_uhl() (matsubara_s_*).
      //=====================================================================

      extern int n_bath_modes;   // Number of auxiliary bath modes per atom

      //---------------------------------------------------------------------
      // Log-bath λ-range search tuning knobs (quantum_no_zero only)
      //
      //   bath_scan_resolution  — grid points per axis for the (log_ls,
      //                           log_le) scan.  Total objective evaluations
      //                           = bath_scan_resolution² · bath_scan_omega_points.
      //   bath_scan_omega_points — points in the linear ω grid used to
      //                           evaluate the SSE objective.  Grid spans
      //                           [ω_max/N, ω_max] with ω_max = 15·max(T, ω₀).
      //   bath_scan_decades     — width (in log-decades) of the asymmetric
      //                           search box around log10(T) and
      //                           log10(max(T, ω₀)).
      //
      // Defaults reproduce the cmp_noise reference algorithm out of the box.
      //---------------------------------------------------------------------
      extern int    bath_scan_resolution;     // default 50
      extern int    bath_scan_omega_points;   // default 500
      extern double bath_scan_decades;        // default 3.5

      extern std::vector<double> matsubara_s_x;
      extern std::vector<double> matsubara_s_y;
      extern std::vector<double> matsubara_s_z;

      // Per-atom noise drawn once per RK4 step (reused across K1-K4)
      extern std::vector<double> qn_x_array;
      extern std::vector<double> qn_y_array;
      extern std::vector<double> qn_z_array;

      // White zero-mode amplitude sqrt(2*T/dt) — used by generate_quantum_noise_orn_uhl.
      // The quantum noise generators use mtrandom::gaussian() (Ziggurat, ~4× faster
      // than Box-Muller) — the same RNG every other Vampire module uses, kept
      // rank-correlated by Vampire's parallel_rng_seed setup at startup.
      extern double mats_white_amp;

      //=====================================================================
      // Exact Ornstein-Uhlenbeck coefficients (quantum_zero noise type)
      //
      //   s_n(t+dt) = ou_decay[n]*s_n(t) + ou_diffuse[n]*xi_n
      //   noise contribution = ou_drift[n]*s_n_old + ou_noise_amp[n]*xi_n
      //=====================================================================
      extern std::vector<double> ou_decay;      // exp(-nu_n * dt)
      extern std::vector<double> ou_diffuse;    // sqrt(1 - decay_n^2)
      extern std::vector<double> ou_drift;      // sqrt(2T/nu_n) * (decay - 1) / dt
      extern std::vector<double> ou_noise_amp;  // sqrt(2T/nu_n) * diffuse / dt

      //=====================================================================
      // Optimised-range log-bath state and coefficients (quantum_no_zero)
      //
      // K = n_bath_modes auxiliary OU processes with log-spaced rates lambda_k.
      // Lambda range [lambda_lo, lambda_hi] found at init by 2-D grid scan
      // (see noise_ho.cpp / cmp_noise.cpp:lb_opt_scan) minimising squared
      // error against omega*(coth(omega/2T)-1) on a linear ω grid.
      // Weight formula: w_k = sqrt((coth(lam_k/2T)-1) * dl_k / lam_k).
      //=====================================================================

      // Per-atom state arrays: size = num_atoms * n_bath_modes
      extern std::vector<double> lb_opt_s_x;
      extern std::vector<double> lb_opt_s_y;
      extern std::vector<double> lb_opt_s_z;
      // λ-range frozen at init (grid-scan output). Held fixed across T changes —
      // re-running the scan per step is too expensive, but the resulting bath
      // remains an acceptable fit for modest T excursions.
      extern std::vector<double> lb_opt_lambdas;  // size n_bath_modes
      // Precomputed per-mode coefficients (refreshed when sim::temperature changes)
      extern std::vector<double> lb_opt_decay;    // exp(-lambda_k * dt)
      extern std::vector<double> lb_opt_coeff_s;  // sqrt(1 - decay_k^2)
      extern std::vector<double> lb_opt_amp;      // renormalized w_k

      //=====================================================================
      // Cascade OU bath (quantum_no_zero only).
      //
      // Each log-bath mode k drives a chain of n_cascade_modes OU processes:
      //   s_{k,1}  <-- white noise (standard OU, same as existing)
      //   s_{k,j}  <-- s_{k,j-1}  (j=2..M), update: decay_k*s + casc_coeff_k*s_prev
      //
      // The output field is built from the last cascade level s_{k,M}.
      // n_cascade_modes = 1 disables the cascade (identical to existing path).
      // Set via quantum:cascade-modes.
      //=====================================================================
      extern int n_cascade_modes;

      // State layout: [n_atoms * n_bath_modes * n_cascade_modes], indexed as
      //   [atom * K * M  +  mode * M  +  level]
      extern std::vector<double> lb_opt_casc_s_x;
      extern std::vector<double> lb_opt_casc_s_y;
      extern std::vector<double> lb_opt_casc_s_z;

      // Per-mode coupling coefficient (1 - decay_k) / sqrt(lambda_k). Fixed at setup.
      extern std::vector<double> lb_opt_casc_coeff;

      void generate_quantum_noise_log_bath_opt_cascade(const int atom, const double amp, const double dt);

      //=====================================================================
      // 4th-order Butterworth low-pass post-filter (quantum_no_zero only).
      //
      // Wraps generate_quantum_noise_log_bath_opt: calls the raw generator
      // first, then passes the result through two cascaded biquad sections
      // to suppress the 1/ω² high-frequency tail without touching bath state.
      //
      // Filter design: bilinear transform with cutoff pre-warping.
      //   4th-order Butterworth  =  stage 0 (Q ≈ 1.3066)  ×  stage 1 (Q ≈ 0.5412)
      //
      // Activation: set quantum:filter-cutoff = <value> !THz in the input file.
      // When butter_cutoff_Hz == 0 the filter path is never entered.
      //=====================================================================

      // Coefficients for one 2nd-order IIR section (a0 normalised to 1).
      struct BiquadCoeffs {
         double b0, b1, b2;   // FIR (numerator) coefficients
         double a1, a2;        // IIR (denominator) coefficients
      };

      // Two biquad sections that together form the 4th-order Butterworth.
      extern BiquadCoeffs butter_stage[2];

      // Becomes true after coefficients are first computed; reset on cutoff change.
      extern bool   butter_coeffs_valid;

      // Cutoff stored by setup_log_bath_filtered; used at first generate call
      // (dt is available there as a parameter, not at setup time).
      extern double butter_cutoff_Hz;

      // Per-atom filter state.  Flat 1D layout: [atom * 4 + stage * 2 + delay].
      //   stage  ∈ {0, 1}  — the two cascaded biquad sections
      //   delay  ∈ {0, 1}  — s1, s2 of the Transposed Direct Form II
      // Allocated in setup_log_bath_filtered; one vector per Cartesian axis.
      extern std::vector<double> filter_state_x;
      extern std::vector<double> filter_state_y;
      extern std::vector<double> filter_state_z;

      // Compute bilinear-transform biquad coefficients for one 2nd-order section.
      //   dt_s       — simulation time step in seconds
      //   cutoff_Hz  — desired -3 dB cutoff frequency in Hz
      //   Q          — quality factor of this Butterworth section
      BiquadCoeffs calculate_biquad_coeffs(double dt_s, double cutoff_Hz, double Q);

      // Initialise the log-bath AND the Butterworth filter state.
      // Must be called instead of (not in addition to) setup_log_bath_opt.
      void setup_log_bath_filtered(int num_atoms_total, double cutoff_frequency_THz);

      // Generate quantum noise and immediately LP-filter it in place.
      // Overwrites qn_x/y/z_array[atom] with the filtered field.
      void generate_quantum_noise_log_bath_opt_filtered(const int atom, const double amp, const double dt);

      //=====================================================================
      // Cached T_scaled used by the last refresh of the quantum-noise
      // coefficient arrays.  refresh_quantum_noise_for_T() compares the
      // current scale_temperature(sim::temperature) against this and only
      // re-runs the refresh helpers when they differ.  Sentinel value
      // -1.0 forces a refresh on the first step.
      //=====================================================================
      extern double last_T_scaled;

      //=====================================================================
      // SLD phonon quantum noise — parallel OU bath for the lattice subsystem.
      //
      // Independent OU oscillator states for the phononic degrees of freedom.
      // Generated once per Suzuki-Trotter step and stored in qn_phonon_*_array.
      // The classical amplitude reference is F_th_sigma (force noise prefactor).
      //=====================================================================

      // Per-atom output noise arrays for the phonon (lattice) subsystem
      extern std::vector<double> qn_phonon_x_array;
      extern std::vector<double> qn_phonon_y_array;
      extern std::vector<double> qn_phonon_z_array;

      // Per-atom OU bath state for phonon noise (Matsubara modes, quantum_zero)
      extern std::vector<double> matsubara_phonon_s_x;
      extern std::vector<double> matsubara_phonon_s_y;
      extern std::vector<double> matsubara_phonon_s_z;

      // Per-atom OU bath state for phonon noise (log-bath modes, quantum_no_zero)
      extern std::vector<double> lb_opt_phonon_s_x;
      extern std::vector<double> lb_opt_phonon_s_y;
      extern std::vector<double> lb_opt_phonon_s_z;

      // Per-material phonon OU coefficients (Matsubara path)
      extern std::vector<double> ou_phonon_decay;
      extern std::vector<double> ou_phonon_diffuse;
      extern std::vector<double> ou_phonon_drift;
      extern std::vector<double> ou_phonon_noise_amp;
      extern double mats_phonon_white_amp;

      // Per-material phonon OU coefficients (log-bath path)
      extern std::vector<double> lb_opt_phonon_decay;
      extern std::vector<double> lb_opt_phonon_coeff_s;
      extern std::vector<double> lb_opt_phonon_amp;
      extern std::vector<double> lb_opt_phonon_lambdas;

      // Per-material phonon amplitude and frequency parameters (computed at init)
      extern std::vector<double> material_phonon_amp_array;    // overall force amplitude
      extern std::vector<double> material_phonon_omega0_array; // characteristic phonon frequency
      extern std::vector<double> material_phonon_gamma_array;  // phonon damping rate

      // Per-material SLD spin-noise amplitude (computed at init from H_th_sigma so the
      // high-T white limit reproduces the classical SLD spin noise H_th_sigma*sqrt(T)).
      extern std::vector<double> material_sld_spin_amp_array;

      // Per-material classical spin-noise sigma = H_th_sigma (no spin_pref factor).
      // Used only by the classical branch of sld_noise::generate().
      extern std::vector<double> material_sld_classical_sigma_array;

      // Cached T for phonon refresh (separate from spin last_T_scaled)
      extern double last_T_scaled_phonon;

      //=====================================================================
      // SLD / ASD FFT spin noise (pre-generated, pure coth / coth-1 shape)
      //
      // Flat layout: [atom * sld_fft_n_fine + step] for each component.
      // Allocated and filled once by generate_sld_spin_fft() called from
      // sld_noise::initialize() when kind is quantum_fft or quantum_no_zero_fft.
      // sld_fft_step_index is advanced by sld_noise::generate() each step.
      //=====================================================================
      extern std::vector<double> sld_fft_spin_x;  // [num_atoms * n_fine]
      extern std::vector<double> sld_fft_spin_y;
      extern std::vector<double> sld_fft_spin_z;
      extern uint64_t sld_fft_step_index;          // current fine step
      extern uint64_t sld_fft_n_fine;              // total steps pre-generated
      extern int      sld_fft_num_atoms;            // atoms used at generation time

      //=====================================================================
      // SLD / ASD spin-noise export state (set by sld_noise::enable_spin_noise_export)
      //=====================================================================
      extern bool        sld_export_noise;
      extern std::string sld_export_noise_filename;
      extern int         sld_export_noise_atom;

      //=====================================================================
      // Internal function declarations
      //=====================================================================

      // --- LLG time stepper functions ---
      void llg_HO();             // Serial HO stepper
      void llg_HO_mpi();         // MPI-parallel HO stepper
      void llg_FFT();            // Serial FFT stepper
      void llg_FFT_mpi();        // MPI-parallel FFT stepper

      // --- Equations of motion ---

      /// HO method: noise_x/y/z are pre-drawn and held fixed across K1-K4
      void LL_HO_method(const double* y, const double* H, double* dydt,
                        const int material, const double dt,
                        const double noise_x, const double noise_y, const double noise_z);

      /// FFT method: same equations but noise enters through H, not the momentum equation
      void LL_FFT_method(const double* y, const double* H, double* dydt,
                         const int material);

      // --- Per-atom RK4 primitives (shared by all four time-steppers) ---

      /// Save current (S, q, p) state into y_in_storage[atom]
      void save_initial_state(const int atom);

      /// Normalise the spin part (first three components) of a 9-vector in place
      void renormalize_spin(double* y);

      /// HO effective field: H = spin_field + external_field
      void collect_H_HO(const int atom, double H[3]);

      /// FFT effective field: H = spin_field + external_field + qn_{x,y,z}_array[atom]
      /// Noise is read from qn_*_array, populated once per step by
      /// draw_noise_all_atoms_FFT(), and reused across K1-K4.
      void collect_H_FFT(const int atom, double H[3]);

      /// y_pred = y_in + coeff*k_stage; renormalise spin part of y_pred
      void predict_and_renorm(const int atom, const double* k_stage, const double coeff);

      /// y_pred = y_in + (dt/6)*(k1+2k2+2k3+k4); renormalise spin part of y_pred
      void rk4_combine_and_renorm(const int atom, const double dt_over_6);

      /// Copy y_pred[0..2] back to atoms::*_spin_array (used between RK4 stages)
      void writeback_spin(const int atom);

      /// Copy all 9 components of y_pred back to atoms + q/p arrays
      void writeback_full_state(const int atom);

      /// Draw thermal noise for atoms[lo..hi) (HO method). Classical: white +
      /// FDT prefactor; quantum: dispatch to generate_quantum_noise_HO().
      void draw_noise_all_atoms_HO(const int lo, const int hi, const double dt);

      /// Sample the FFT noise field for atoms[lo..hi) into qn_{x,y,z}_array.
      /// One sample per atom per RK4 step, held fixed across K1-K4.
      void draw_noise_all_atoms_FFT(const int lo, const int hi);

      // --- Quantum noise generation ---

      /// Dispatch to the correct generator based on noise_type; stores result in qn_x/y/z_array[atom]
      void generate_quantum_noise_HO(const int atom, const int material, const double dt);

      /// Exact Ornstein-Uhlenbeck (unconditionally stable) — quantum_zero noise type.
      /// amp is the overall noise amplitude (thermostat: sqrt(gamma*A/S0); SLD spin:
      /// material_sld_spin_amp_array), kept as a parameter so the same bath machinery
      /// serves both the thermostat and the SLD spin subsystem.
      void generate_quantum_noise_orn_uhl(const int atom, const double amp, const double dt);

      /// Log-spaced OU bath, lambda range optimised by 2-D grid scan — quantum_no_zero noise type
      void generate_quantum_noise_log_bath_opt(const int atom, const double amp, const double dt);

      /// Squared-error objective for the log-bath lambda-range grid scan (shared by the
      /// thermostat spin setup and the SLD phonon/spin setup).
      double lb_objective(const double log_ls, const double log_le,
                          const double T_scaled, const int K,
                          const std::vector<double>& scan_om);

      // --- Quantum noise setup (called from initialize_HO) ---

      /// Allocate Matsubara state arrays and precompute OU coefficients (quantum_zero)
      void setup_orn_uhl(const int num_atoms_total);

      /// Allocate log-bath state arrays, run the 2-D grid scan to choose the
      /// optimal lambda range, precompute log-bath coefficients (quantum_no_zero)
      void setup_log_bath_opt(const int num_atoms_total);

      // --- Dynamic-T refresh (called at the top of each LLG step) ---

      /// Recompute the T-dependent OU coefficient arrays at the given T.
      /// Cheap: O(n_bath_modes).
      void refresh_orn_uhl_coeffs(const double T_scaled, const double dt);

      /// Recompute the T-dependent log-bath amplitude arrays at the given T.
      /// Uses the fixed λ-range from the init-time grid scan.
      /// Cheap: O(n_bath_modes).
      void refresh_log_bath_opt_coeffs(const double T_scaled, const double dt);

      /// Refresh whichever quantum-noise coefficients are active *only if*
      /// scale_temperature(sim::temperature) differs from last_T_scaled.
      /// Called at the top of llg_HO / llg_HO_mpi; no-op when T is constant.
      void refresh_quantum_noise_for_T();

      // --- SLD quantum noise (defined in noise_sld.cpp; driven via quantum::sld_noise) ---

      /// Generate phonon quantum noise for one atom; stores result in qn_phonon_*_array[atom].
      void generate_quantum_noise_phonon_HO(const int atom, const int material, const double dt);

      /// Refresh phonon OU coefficients if temperature changed since last call.
      void refresh_quantum_noise_phonon_for_T();

      /// Generate SLD spin quantum noise for one atom; reuses the thermostat spin bath
      /// machinery with the SLD-derived amplitude (material_sld_spin_amp_array) so the
      /// high-T white limit reproduces the classical SLD spin noise H_th_sigma*sqrt(T).
      /// Result stored in qn_x/y/z_array[atom].
      void generate_sld_spin_HO(const int atom, const int material, const double dt);

      // --- Noise generation and retrieval ---

      /// Compute the power spectral density at frequency omega for temperature T
      double PSD(const double omega, const double T, const int material);

      /// Map atom indices to offsets in the flat coarse_noise_field array
      void assign_unique_indices(int n_coarse, int num_atoms_local);

      /// Shared FFT shaping body used by calculate_noise() and generate_noise_window():
      ///   forward FFT on in_buf → multiply each frequency bin by sqrt_psd[i] → inverse FFT into out_time.
      /// Caller fills in_buf and post-processes out_time (scaling, valid-region extraction).
      void run_fft_pipeline(double* in_buf,
                            fftw_plan plan_fwd, fftw_complex* out_freq,
                            fftw_plan plan_bwd, double* out_time,
                            const std::vector<double>& sqrt_psd, int nc);

      /// Non-windowed: generate all noise at once using a temporary FFT
      void calculate_noise(int realizations, double dt_fine,
                           int M, double T, int n_coarse_total,
                           std::vector<double>& noise_field);

      /// Linearly interpolate a coarse noise value at a fractional fine step index
      double get_noise(const std::vector<double>& coarse_noise,
                       double fine_step_idx, int M, size_t atom_idx);

      /// Export / append actual simulation noise for atom 0 to file (FFT method)
      void export_noise_data(const std::string& filename);

      /// Append one row to the HO noise file: (t, q_x[0], q_y[0], q_z[0]).
      /// q is the auxiliary oscillator that couples the noise to the spin
      /// (dS/dt = S × (H + q)), so q[0](t) is the HO analogue of what
      /// export_noise_data writes for the FFT path.
      void export_ho_noise_step();

      // --- Windowed noise management ---

      /// Allocate persistent FFT resources and precompute scaling constants
      void init_noise_structures(int n_coarse, int realizations,
                                 double dt_fine, int M, double T);

      /// Generate one window of noise using the persistent FFT resources
      void generate_noise_window();

      /// Check if noise_index has exceeded the window; regenerate if so
      void update_noise_if_needed();

      /// Free persistent FFT resources (called at program exit)
      void cleanup_noise_structures();

   } // end of internal namespace

} // end of quantum namespace

#endif //QUANTUM_INTERNAL_H_
