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
//   Quantum noise for the spin-lattice dynamics (SLD) module.
//
//   This file lives in the quantum module so that the colored-noise bath
//   machinery is defined exactly once (the thermostat spin generators live in
//   noise_ho.cpp; this file adds the SLD phonon generators and reuses the spin
//   generators with an SLD-derived amplitude). The SLD integrator talks to it
//   only through the public quantum::sld_noise API below — it never touches
//   quantum::internal.
//
//   Two coupled OU baths are produced per Suzuki-Trotter step:
//
//     * spin   -> qn_{x,y,z}_array        (added to the effective field H_eff)
//                 reuses the thermostat spin bath (Matsubara / log-bath) with
//                 amplitude material_sld_spin_amp_array, derived so the high-T
//                 white limit reproduces SLD's classical spin noise
//                 H_th_sigma * sqrt(T).
//
//     * phonon -> qn_phonon_{x,y,z}_array (added to the atomic force/velocity)
//                 independent OU bath on the phonon time scale (dt_ps), with
//                 amplitude material_phonon_amp_array, derived so the high-T
//                 white limit reproduces F_th_sigma * sqrt(T).
//
//   Only the universal coth(w/2T) (quantum) / coth(w/2T)-1 (quantum-no-zero)
//   spectral *shape* is borrowed from the quantum module; all amplitudes come
//   from the SLD material parameters passed in via sld_noise::initialize().
//   No quantum-lorentzian-* inputs are required for an SLD run.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#ifdef FFT
#include <fftw3.h>
#endif

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"

// Module headers
#include "internal.hpp"

namespace quantum {
namespace internal {

// Toggle: true = batch-fill noise buffer via std::generate before OU loop; false = generate inline
static bool use_preallocated_noise = true;

   //==========================================================================
   // SLD spin quantum noise
   //
   // Reuses the thermostat spin bath state and coefficients (matsubara_s_* /
   // lb_opt_s_* + ou_* / lb_opt_*), which are allocated by setup_orn_uhl /
   // setup_log_bath_opt and refreshed by refresh_quantum_noise_for_T (both on
   // the spin time scale, scale_temperature + mp::dt). Only the amplitude
   // differs from the thermostat path: material_sld_spin_amp_array instead of
   // sqrt(gamma*A/S0). Result stored in qn_{x,y,z}_array.
   //==========================================================================
   void generate_sld_spin_HO(const int atom, const int material, const double dt) {
      const double amp = material_sld_spin_amp_array[material];
      switch (noise_type) {
         case quantum_zero:    generate_quantum_noise_orn_uhl(atom, amp, dt);        break;
         case quantum_no_zero:
            if (n_cascade_modes > 1)
               generate_quantum_noise_log_bath_opt_cascade(atom, amp, dt);
            else
               generate_quantum_noise_log_bath_opt(atom, amp, dt);
            break;
         default: break;
      }
   }

   //==========================================================================
   // SLD phonon quantum noise — parallel OU bath for the lattice subsystem.
   //
   // Independent per-atom state arrays (matsubara_phonon_s_* / lb_opt_phonon_s_*),
   // output stored in qn_phonon_{x,y,z}_array, on the phonon time scale
   // dt_ps = mp::dt_SI*1e12 and T_scaled = scale_temperature_phonon(T) (rad/ps).
   // Amplitude material_phonon_amp_array matches F_th_sigma*sqrt(T) at high T.
   //==========================================================================

   void refresh_orn_uhl_phonon_coeffs(const double T_scaled, const double dt) {
      const double inv_dt = 1.0 / dt;
      for (int i = 0; i < n_bath_modes; ++i) {
         const double nu        = 2.0 * M_PI * (i + 1) * T_scaled;
         const double decay_i   = std::exp(-nu * dt);
         const double diffuse_i = std::sqrt(1.0 - decay_i * decay_i);
         ou_phonon_decay[i]     = decay_i;
         ou_phonon_diffuse[i]   = diffuse_i;
         ou_phonon_drift[i]     = std::sqrt(2.0 * T_scaled / nu) * (decay_i - 1.0) * inv_dt;
         ou_phonon_noise_amp[i] = std::sqrt(2.0 * T_scaled / nu) * diffuse_i * inv_dt;
      }
      mats_phonon_white_amp = std::sqrt(2.0 * T_scaled * inv_dt);
   }

   void refresh_log_bath_opt_phonon_coeffs(const double T_scaled, const double dt) {
      const int K = n_bath_modes;
      const double* lam = lb_opt_phonon_lambdas.data();

      std::vector<double> raw_amp(K, 0.0);
      for (int k = 0; k < K - 1; ++k) {
         const double dl = lam[k + 1] - lam[k];
         const double x  = lam[k] / (2.0 * T_scaled);
         const double ct = (std::fabs(x) < 1e-10) ? 1.0 / x : 1.0 / std::tanh(x);
         raw_amp[k] = std::sqrt(lam[k] * std::max(0.0, ct - 1.0) * dl);
      }

      const double omega_ref = material_phonon_omega0_array.empty()
                               ? T_scaled
                               : material_phonon_omega0_array[0];
      double S_model = 0.0;
      for (int k = 0; k < K; ++k)
         S_model += 2.0 * lam[k] * raw_amp[k] * raw_amp[k] / (omega_ref * omega_ref + lam[k] * lam[k]);
      const double x_ref = omega_ref / (2.0 * T_scaled);
      const double ct_ref = (std::fabs(x_ref) < 1e-10) ? 1.0 / x_ref : 1.0 / std::tanh(x_ref);
      const double S_tgt  = omega_ref * std::max(0.0, ct_ref - 1.0);
      const double renorm = (S_model > 0.0) ? std::sqrt(S_tgt / S_model) : 1.0;

      for (int k = 0; k < K; ++k) {
         const double d = std::exp(-lam[k] * dt);
         lb_opt_phonon_decay[k]   = d;
         lb_opt_phonon_coeff_s[k] = std::sqrt(1.0 - d * d);
         lb_opt_phonon_amp[k]     = raw_amp[k] * renorm;
      }
   }

   void generate_quantum_noise_phonon_orn_uhl(const int atom, const int material, const double /*dt*/) {
      const double amp  = material_phonon_amp_array[material];
      const int    nm   = n_bath_modes;
      const size_t base = static_cast<size_t>(atom) * nm;

      double aux_x = 0.0, aux_y = 0.0, aux_z = 0.0;

      thread_local static std::vector<double> noise_buf_ph;
      if (use_preallocated_noise) {
         noise_buf_ph.resize(static_cast<size_t>(nm) * 3);
         std::generate(noise_buf_ph.begin(), noise_buf_ph.end(), mtrandom::gaussian);
      }

      int i = 0;
      for (; i + 1 < nm; i += 2) {
         const double dc0 = ou_phonon_decay[i],     dc1 = ou_phonon_decay[i + 1];
         const double df0 = ou_phonon_diffuse[i],   df1 = ou_phonon_diffuse[i + 1];
         const double dr0 = ou_phonon_drift[i],     dr1 = ou_phonon_drift[i + 1];
         const double na0 = ou_phonon_noise_amp[i], na1 = ou_phonon_noise_amp[i + 1];

         double& sx0 = matsubara_phonon_s_x[base + i];  double& sx1 = matsubara_phonon_s_x[base + i + 1];
         double& sy0 = matsubara_phonon_s_y[base + i];  double& sy1 = matsubara_phonon_s_y[base + i + 1];
         double& sz0 = matsubara_phonon_s_z[base + i];  double& sz1 = matsubara_phonon_s_z[base + i + 1];

         double zx0, zx1, zy0, zy1, zz0, zz1;
         if (use_preallocated_noise) {
            zx0 = noise_buf_ph[i];          zx1 = noise_buf_ph[i + 1];
            zy0 = noise_buf_ph[nm + i];     zy1 = noise_buf_ph[nm + i + 1];
            zz0 = noise_buf_ph[2 * nm + i]; zz1 = noise_buf_ph[2 * nm + i + 1];
         } else {
            zx0 = mtrandom::gaussian();  zx1 = mtrandom::gaussian();
            zy0 = mtrandom::gaussian();  zy1 = mtrandom::gaussian();
            zz0 = mtrandom::gaussian();  zz1 = mtrandom::gaussian();
         }

         aux_x += dr0 * sx0 + na0 * zx0 + dr1 * sx1 + na1 * zx1;
         aux_y += dr0 * sy0 + na0 * zy0 + dr1 * sy1 + na1 * zy1;
         aux_z += dr0 * sz0 + na0 * zz0 + dr1 * sz1 + na1 * zz1;

         sx0 = dc0 * sx0 + df0 * zx0;  sx1 = dc1 * sx1 + df1 * zx1;
         sy0 = dc0 * sy0 + df0 * zy0;  sy1 = dc1 * sy1 + df1 * zy1;
         sz0 = dc0 * sz0 + df0 * zz0;  sz1 = dc1 * sz1 + df1 * zz1;
      }
      if (i < nm) {
         const double dc = ou_phonon_decay[i],  df = ou_phonon_diffuse[i];
         const double dr = ou_phonon_drift[i],  na = ou_phonon_noise_amp[i];
         double zx, zy, zz;
         if (use_preallocated_noise) {
            zx = noise_buf_ph[i]; zy = noise_buf_ph[nm + i]; zz = noise_buf_ph[2 * nm + i];
         } else { zx = mtrandom::gaussian(); zy = mtrandom::gaussian(); zz = mtrandom::gaussian(); }
         aux_x += dr * matsubara_phonon_s_x[base + i] + na * zx;
         aux_y += dr * matsubara_phonon_s_y[base + i] + na * zy;
         aux_z += dr * matsubara_phonon_s_z[base + i] + na * zz;
         matsubara_phonon_s_x[base + i] = dc * matsubara_phonon_s_x[base + i] + df * zx;
         matsubara_phonon_s_y[base + i] = dc * matsubara_phonon_s_y[base + i] + df * zy;
         matsubara_phonon_s_z[base + i] = dc * matsubara_phonon_s_z[base + i] + df * zz;
      }

      const double w_x = mtrandom::gaussian();
      const double w_y = mtrandom::gaussian();
      const double w_z = mtrandom::gaussian();

      qn_phonon_x_array[atom] = amp * (mats_phonon_white_amp * w_x + aux_x);
      qn_phonon_y_array[atom] = amp * (mats_phonon_white_amp * w_y + aux_y);
      qn_phonon_z_array[atom] = amp * (mats_phonon_white_amp * w_z + aux_z);
   }

   void generate_quantum_noise_phonon_log_bath_opt(const int atom, const int material, const double /*dt*/) {
      const double amp  = material_phonon_amp_array[material];
      const int    K    = n_bath_modes;
      const size_t base = static_cast<size_t>(atom) * K;

      double aux_x = 0.0, aux_y = 0.0, aux_z = 0.0;
      for (int k = 0; k < K; ++k) {
         const double d  = lb_opt_phonon_decay[k];
         const double cs = lb_opt_phonon_coeff_s[k];
         const double a  = lb_opt_phonon_amp[k];
         const double zx = mtrandom::gaussian();
         const double zy = mtrandom::gaussian();
         const double zz = mtrandom::gaussian();
         lb_opt_phonon_s_x[base + k] = d * lb_opt_phonon_s_x[base + k] + cs * zx;
         lb_opt_phonon_s_y[base + k] = d * lb_opt_phonon_s_y[base + k] + cs * zy;
         lb_opt_phonon_s_z[base + k] = d * lb_opt_phonon_s_z[base + k] + cs * zz;
         aux_x += a * lb_opt_phonon_s_x[base + k];
         aux_y += a * lb_opt_phonon_s_y[base + k];
         aux_z += a * lb_opt_phonon_s_z[base + k];
      }
      qn_phonon_x_array[atom] = amp * aux_x;
      qn_phonon_y_array[atom] = amp * aux_y;
      qn_phonon_z_array[atom] = amp * aux_z;
   }

   void generate_quantum_noise_phonon_HO(const int atom, const int material, const double dt) {
      switch (noise_type) {
         case quantum_zero:    generate_quantum_noise_phonon_orn_uhl(atom, material, dt);       break;
         case quantum_no_zero: generate_quantum_noise_phonon_log_bath_opt(atom, material, dt);  break;
         default: break;
      }
   }

   void refresh_quantum_noise_phonon_for_T() {
      const double T_now = scale_temperature_phonon(sim::temperature);
      if (T_now == last_T_scaled_phonon) return;
      const double dt = mp::dt_SI * 1e12;  // phonon time step in ps
      if (noise_type == quantum_zero) {
         refresh_orn_uhl_phonon_coeffs(T_now, dt);
      } else if (noise_type == quantum_no_zero) {
         refresh_log_bath_opt_phonon_coeffs(T_now, dt);
      }
      last_T_scaled_phonon = T_now;
   }

   //==========================================================================
   // FFT-based SLD/ASD spin noise pre-generation (pure coth / coth-1 shape)
   //
   // Called once from sld_noise::initialize() when kind == quantum_fft or
   // quantum_no_zero_fft.  Fills sld_fft_spin_{x,y,z} with coloured noise
   // whose spectral density is:
   //   quantum_fft     : PSD(k) ~ H_th_sigma^2 * coth(omega_k / 2T_scaled)
   //   quantum_no_zero_fft: PSD(k) ~ H_th_sigma^2 * (coth(omega_k/2T_scaled) - 1)
   //
   // Amplitude calibration: the same H_th_sigma * spin_pref factor as the
   // HO path, so both paths agree on the per-frequency amplitude.
   //
   // Storage layout: sld_fft_spin_x[ atom * n_fine + step ]
   //==========================================================================
   static void generate_sld_spin_fft(int num_atoms, uint64_t n_fine,
                                      bool no_zero,
                                      const std::vector<double>& H_th_sigma_per_mat) {

#ifdef FFT
      using namespace quantum::internal;

      const double dt       = mp::dt;
      const double T_scaled = scale_temperature(sim::temperature);
      const double spin_pref = std::sqrt(dt * constants::hbar
                                         * constants::gyromagnetic_ratio
                                         / (2.0 * constants::kB));

      const size_t total = static_cast<size_t>(num_atoms) * n_fine;

      // Memory guard: warn and abort FFT if > 8 GiB
      const size_t bytes = total * 3 * sizeof(double);
      if (bytes > 8UL * 1024 * 1024 * 1024) {
         std::cerr << "Warning: FFT spin-noise array would require "
                   << bytes / (1024*1024) << " MiB — exceeds 8 GiB limit.\n"
                   << "  Falling back to HO (on-the-fly) noise generation.\n"
                   << "  Reduce system size or total steps, or use noise-type=quantum-no-zero.\n";
         sld_fft_n_fine    = 0;
         sld_fft_num_atoms = 0;
         return;
      }

      sld_fft_spin_x.assign(total, 0.0);
      sld_fft_spin_y.assign(total, 0.0);
      sld_fft_spin_z.assign(total, 0.0);
      sld_fft_n_fine    = n_fine;
      sld_fft_num_atoms = num_atoms;
      sld_fft_step_index = 0;

      const int nc = static_cast<int>(n_fine);

      // One FFTW plan pair (r2c + c2r) reused for every (atom, component)
      double*        in  = static_cast<double*>(fftw_malloc(sizeof(double) * nc));
      fftw_complex*  out = static_cast<fftw_complex*>(
                              fftw_malloc(sizeof(fftw_complex) * (nc/2 + 1)));
      fftw_plan fwd = fftw_plan_dft_r2c_1d(nc, in,  out, FFTW_ESTIMATE);
      fftw_plan bwd = fftw_plan_dft_c2r_1d(nc, out, in,  FFTW_ESTIMATE);

      // sqrt of the target PSD, DC bin = 0 (no mean shift).
      //
      // The targets are  omega*coth(omega/2T)  and  omega*(coth(omega/2T)-1),
      // NOT coth and coth-1: the leading factor of omega is part of the
      // spectrum, not a normalisation. Without it the generated noise falls as
      // 1/omega relative to the target, which at 300 K leaves only ~3% of the
      // correct power over 0.1-3 T_scaled and breaks the classical DC limit
      // that both spectra must satisfy.
      //
      // With the factor restored both targets tend to the classical plateau
      // 2*T_scaled as omega -> 0, matching the white-noise amplitude exactly,
      // and neither diverges at the origin.
      //
      // coth(x)-1 is evaluated as 2/expm1(2x): the literal difference loses all
      // precision for x > ~18, where the true value is ~1e-16.
      std::vector<double> sqrt_psd(nc/2 + 1, 0.0);
      for (int k = 1; k <= nc/2; ++k) {
         const double omega = 2.0 * M_PI * k / (nc * dt);
         const double x     = omega / (2.0 * T_scaled);
         double psd;
         if (no_zero) {
            // omega * (coth(x) - 1) = 2*omega / (exp(2x) - 1)
            psd = (x > 350.0) ? 0.0 : 2.0 * omega / std::expm1(2.0 * x);
         } else {
            // omega * coth(x), with the omega/x -> 2T limit as x -> 0
            const double coth = (std::fabs(x) < 1e-10) ? 1.0 / x : 1.0 / std::tanh(x);
            psd = omega * coth;
         }
         sqrt_psd[k] = std::sqrt(std::max(0.0, psd));
      }

      // white-noise sigma = 1/sqrt(dt); final_scale = amp / n_fine * sigma
      const double sigma_white = 1.0 / std::sqrt(dt);
      const double inv_nc      = 1.0 / nc;

      std::cout << "  Pre-generating FFT spin noise for " << num_atoms
                << " atoms × " << n_fine << " steps..." << std::endl;

      for (int atom = 0; atom < num_atoms; ++atom) {
         const int    mat = atoms::type_array[atom];
         const double amp = H_th_sigma_per_mat[mat] * spin_pref;
         const double fsc = amp * sigma_white * inv_nc;   // per-sample final scale

         const size_t base = static_cast<size_t>(atom) * n_fine;

         auto run_component = [&](std::vector<double>& dst) {
            for (int n = 0; n < nc; ++n) in[n] = mtrandom::gaussian();
            fftw_execute(fwd);
            for (int k = 0; k <= nc/2; ++k) {
               out[k][0] *= sqrt_psd[k];
               out[k][1] *= sqrt_psd[k];
            }
            fftw_execute(bwd);
            for (int n = 0; n < nc; ++n)
               dst[base + n] = in[n] * fsc;
         };

         run_component(sld_fft_spin_x);
         run_component(sld_fft_spin_y);
         run_component(sld_fft_spin_z);
      }
      std::cout << "  FFT spin-noise pre-generation done." << std::endl;

      fftw_destroy_plan(fwd);
      fftw_destroy_plan(bwd);
      fftw_free(in);
      fftw_free(out);

#else
      std::cerr << "Warning: FFTW not compiled in — quantum_fft noise kind is unavailable.\n"
                << "  Recompile with -DFFT -lfftw3, or use noise-type=quantum / quantum-no-zero.\n";
      sld_fft_n_fine    = 0;
      sld_fft_num_atoms = 0;
      (void)num_atoms; (void)n_fine; (void)no_zero; (void)H_th_sigma_per_mat;
#endif
   }

} // end of internal namespace

   //==========================================================================
   // Public SLD-noise API (quantum::sld_noise) — the only surface the SLD
   // integrator uses. Keeps src/spinlattice free of quantum::internal.
   //==========================================================================
   namespace sld_noise {

      //----------------------------------------------------------------------
      // One-time setup. Maps the SLD noise kind onto the quantum noise_type,
      // derives the spin and phonon amplitudes from the SLD material params,
      // and allocates / scans both baths. Self-contained: no dependency on
      // quantum::initialize() or on any quantum-lorentzian-* input.
      //----------------------------------------------------------------------
      void initialize(kind_t kind, int num_atoms,
                      const std::vector<material_params>& mats,
                      uint64_t n_fine) {

         using namespace quantum::internal;

         // Detect FFT variants
         const bool use_fft = (kind == kind_t::quantum_fft || kind == kind_t::quantum_no_zero_fft);

         // The pre-generated trace is shaped once, at the temperature holding
         // when this runs, and cannot follow a temperature that moves later.
         // The thermostat path (llg-quantum) refuses such programs in
         // supported_program(); this path had no equivalent check, so a laser
         // pulse combined with an -fft noise type ran the whole simulation on
         // noise for the starting temperature without saying so.
         if (use_fft && dynamic_temperature_program()) {
            std::cerr << "Error: pre-generated (FFT) quantum noise cannot be used with a "
                      << "program that varies the temperature.\n"
                      << "  Its spectrum is fixed at the temperature holding when the noise "
                      << "is generated,\n  so the run would silently use the wrong "
                      << "temperature throughout.\n"
                      << "  Use the on-the-fly generator instead: drop the '-fft' suffix "
                      << "from the noise type\n  (quantum-fft -> quantum, "
                      << "quantum-no-zero-fft -> quantum-no-zero)." << std::endl;
            err::vexit();
         }
         const bool no_zero = (kind == kind_t::quantum_no_zero || kind == kind_t::quantum_no_zero_fft);

         // Map SLD kind -> internal thermostat noise_type
         switch (kind) {
            case kind_t::classical:           noise_type = quantum::internal::classical;      break;
            case kind_t::quantum:             noise_type = quantum::internal::quantum_zero;    break;
            case kind_t::quantum_no_zero:     noise_type = quantum::internal::quantum_no_zero; break;
            case kind_t::quantum_fft:         noise_type = quantum::internal::quantum_zero;    break;
            case kind_t::quantum_no_zero_fft: noise_type = quantum::internal::quantum_no_zero; break;
         }

         // Resize spin and phonon output arrays
         qn_x_array.assign(num_atoms, 0.0);
         qn_y_array.assign(num_atoms, 0.0);
         qn_z_array.assign(num_atoms, 0.0);
         qn_phonon_x_array.assign(num_atoms, 0.0);
         qn_phonon_y_array.assign(num_atoms, 0.0);
         qn_phonon_z_array.assign(num_atoms, 0.0);

         const int n_mats = static_cast<int>(mats.size());
         material_phonon_amp_array.resize(n_mats);
         material_phonon_omega0_array.resize(n_mats);
         material_phonon_gamma_array.resize(n_mats);
         material_sld_spin_amp_array.resize(n_mats);
         material_sld_classical_sigma_array.resize(n_mats);

         // hbar in eV*ps — used for phonon amplitude and temperature scaling
         static const double hbar_eVps = constants::hbar * 1e12 / 1.6021766208e-19;

         // Spin amplitude prefactor: amp_spin = H_th_sigma * sqrt(mp::dt * hbar * gamma / (2 kB)).
         // Derivation: spin white-mode = amp_spin * sqrt(2*T_scaled/mp::dt) with
         // T_scaled = T*kB/(hbar*gamma); requiring this to equal the classical
         // spin noise H_th_sigma*sqrt(T) cancels T and leaves a per-material const.
         const double spin_pref = std::sqrt(mp::dt * constants::hbar
                                            * constants::gyromagnetic_ratio
                                            / (2.0 * constants::kB));

         for (int mat = 0; mat < n_mats; ++mat) {
            const double damp_lat   = mats[mat].damp_lat;
            const double mass       = mats[mat].mass;        // [kg]
            const double V0_eV      = mats[mat].V0;          // [eV/Ang^2]
            const double H_th_sigma = mats[mat].H_th_sigma;

            // Phonon amplitude: matches F_th_sigma*sqrt(T) at high T.
            material_phonon_amp_array[mat] = (mass > 0.0)
               ? std::sqrt(damp_lat * hbar_eVps / mass)
               : 0.0;

            // Characteristic phonon frequency from harmonic spring constant V0 [eV/Ang^2]:
            //   omega_ph [rad/ps] = sqrt(V0 * 16.02 / mass) * 1e-12
            material_phonon_omega0_array[mat] = (mass > 0.0 && V0_eV > 0.0)
               ? std::sqrt(V0_eV * 16.02 / mass) * 1e-12
               : scale_temperature_phonon(300.0);

            material_phonon_gamma_array[mat] = damp_lat;

            // Spin amplitude: matches H_th_sigma*sqrt(T) at high T.
            material_sld_spin_amp_array[mat] = H_th_sigma * spin_pref;
            // Classical spin sigma: used directly (no spin_pref factor).
            material_sld_classical_sigma_array[mat] = H_th_sigma;
         }

         const double dt_ph    = mp::dt_SI * 1e12;
         const double T_now    = scale_temperature_phonon(sim::temperature);
         const size_t st_size  = static_cast<size_t>(num_atoms) * n_bath_modes;

         // --- Spin bath setup ---
         // HO kinds: allocate the on-the-fly OU bath now.
         // FFT kinds: deferred until after the pre-generation attempt below,
         // because that attempt may fall back to HO (see there).
         if (!use_fft) {
             if (noise_type == quantum_zero)
               setup_orn_uhl(num_atoms);
            else
               setup_log_bath_opt(num_atoms);
         }

         // --- Phonon bath (always HO, regardless of spin-noise method) ---
         if (noise_type == quantum_zero) {
            matsubara_phonon_s_x.assign(st_size, 0.0);
            matsubara_phonon_s_y.assign(st_size, 0.0);
            matsubara_phonon_s_z.assign(st_size, 0.0);
            ou_phonon_decay.resize(n_bath_modes);
            ou_phonon_diffuse.resize(n_bath_modes);
            ou_phonon_drift.resize(n_bath_modes);
            ou_phonon_noise_amp.resize(n_bath_modes);
            refresh_orn_uhl_phonon_coeffs(T_now, dt_ph);

         } else if (noise_type == quantum_no_zero) {
            lb_opt_phonon_s_x.assign(st_size, 0.0);
            lb_opt_phonon_s_y.assign(st_size, 0.0);
            lb_opt_phonon_s_z.assign(st_size, 0.0);
            lb_opt_phonon_lambdas.resize(n_bath_modes);
            lb_opt_phonon_decay.resize(n_bath_modes);
            lb_opt_phonon_coeff_s.resize(n_bath_modes);
            lb_opt_phonon_amp.resize(n_bath_modes);

            const double omega0_ph = material_phonon_omega0_array.empty()
                                     ? T_now : material_phonon_omega0_array[0];
            const double D      = bath_scan_decades;
            const double ls_min = std::log10(T_now) - D;
            const double ls_max = std::log10(T_now) + 1.0;
            const double le_min = std::log10(std::max(T_now, omega0_ph)) - 0.5;
            const double le_max = std::log10(std::max(T_now, omega0_ph)) + D;

            const int    N_OM   = bath_scan_omega_points;
            const double om_max = std::max(T_now, omega0_ph) * 15.0;
            std::vector<double> scan_om(N_OM);
            for (int ii = 0; ii < N_OM; ++ii) scan_om[ii] = (ii + 1) * om_max / N_OM;

            const int n_scan = bath_scan_resolution;
            double best_err = std::numeric_limits<double>::infinity();
            double best_ls  = ls_min, best_le = le_min;
            for (int i_s = 0; i_s < n_scan; ++i_s) {
               const double log_ls = ls_min + i_s * (ls_max - ls_min) / (n_scan - 1);
               for (int j_s = 0; j_s < n_scan; ++j_s) {
                  const double log_le = le_min + j_s * (le_max - le_min) / (n_scan - 1);
                  if (log_le <= log_ls + 0.05) continue;
                  const double err = lb_objective(log_ls, log_le, T_now, n_bath_modes, scan_om);
                  if (err < best_err) { best_err = err; best_ls = log_ls; best_le = log_le; }
               }
            }
            for (int k = 0; k < n_bath_modes; ++k)
               lb_opt_phonon_lambdas[k] = std::pow(10.0,
                  best_ls + static_cast<double>(k) / (n_bath_modes - 1) * (best_le - best_ls));

            refresh_log_bath_opt_phonon_coeffs(T_now, dt_ph);
         }

         last_T_scaled_phonon = T_now;

         // --- FFT spin noise pre-generation ---
         if (use_fft) {
            if (n_fine == 0) {
               std::cerr << "Warning: quantum FFT noise requested but the total step "
                         << "count is 0; falling back to on-the-fly generation.\n";
            } else {
               std::vector<double> H_th_sigma_vec(static_cast<int>(mats.size()));
               for (int m = 0; m < static_cast<int>(mats.size()); ++m)
                  H_th_sigma_vec[m] = mats[m].H_th_sigma;
               generate_sld_spin_fft(num_atoms, n_fine, no_zero, H_th_sigma_vec);
            }

            // generate_sld_spin_fft() leaves sld_fft_n_fine == 0 when it could
            // not allocate (8 GiB guard) or when FFTW is not compiled in, and
            // announces a fallback to on-the-fly generation. That fallback used
            // to be a lie: the OU bath was never allocated, so the next
            // generate() wrote coefficients into zero-length vectors. Allocate
            // it here so the fallback is real.
            if (sld_fft_n_fine == 0) {
               if (noise_type == quantum_zero) setup_orn_uhl(num_atoms);
               else                            setup_log_bath_opt(num_atoms);
            }
         }

         const char* spin_method = use_fft ? "FFT (pre-generated)" : "HO (on-the-fly)";
         std::cout << "  SLD quantum noise : enabled  (spin " << spin_method
                   << " + phonon HO, " << n_bath_modes << " modes)" << std::endl;
      }

      //----------------------------------------------------------------------
      // Once per Suzuki-Trotter step: refresh T-dependent coefficients then
      // draw one spin + one phonon noise sample per atom (held fixed across
      // all spin/velocity sub-updates of the step).
      //----------------------------------------------------------------------
      void generate(int num_atoms) {
         using namespace quantum::internal;
         if (noise_type == quantum::internal::classical) {
            // Classical white Gaussian noise: amplitude = H_th_sigma * sqrt(T).
            const double sqrt_T = std::sqrt(sim::temperature);
            for (int i = 0; i < num_atoms; ++i) {
               const double amp = material_sld_classical_sigma_array[atoms::type_array[i]] * sqrt_T;
               qn_x_array[i] = amp * mtrandom::gaussian();
               qn_y_array[i] = amp * mtrandom::gaussian();
               qn_z_array[i] = amp * mtrandom::gaussian();
            }
            return;
         }

         // FFT kinds: advance the step pointer; all spin noise already pre-baked.
         // Phonon noise is always HO — generate it on the fly as usual.
         if (sld_fft_n_fine > 0) {
            // Advance only from the second call onward, so that step 0 of the
            // simulation reads pre-generated sample 0 rather than sample 1.
            if (sld_fft_started) ++sld_fft_step_index;
            else                 sld_fft_started = true;

            // Running past the pre-generated trace used to fall through to the
            // never-written qn_*_array, i.e. silently inject ZERO thermal noise
            // for the remainder of the run. Fail loudly instead: the trace is
            // sized from sim::equilibration_time + sim::total_time, so an
            // overrun means the program takes more steps than that.
            if (sld_fft_step_index >= sld_fft_n_fine) {
               std::cerr << "Error: quantum FFT noise exhausted after "
                         << sld_fft_n_fine << " steps.\n"
                         << "  The pre-generated trace is sized from "
                         << "sim:equilibration-time-steps + sim:total-time-steps.\n"
                         << "  Either increase those to cover the whole run, or use the\n"
                         << "  on-the-fly generator (quantum:heun-noise-type = quantum "
                         << "or quantum-no-zero)." << std::endl;
               err::vexit();
            }

            refresh_quantum_noise_phonon_for_T();
            const double dt_ps = mp::dt_SI * 1e12;
            for (int i = 0; i < num_atoms; ++i)
               generate_quantum_noise_phonon_HO(i, atoms::type_array[i], dt_ps);
            return;
         }

         // HO kinds: draw noise for all atoms on the fly.
         refresh_quantum_noise_for_T();
         refresh_quantum_noise_phonon_for_T();

         const double dt_ps = mp::dt_SI * 1e12;
         for (int i = 0; i < num_atoms; ++i) {
            const int mat = atoms::type_array[i];
            generate_sld_spin_HO(i, mat, mp::dt);
            generate_quantum_noise_phonon_HO(i, mat, dt_ps);
         }
      }

      //----------------------------------------------------------------------
      // Per-atom accessors (component 0=x, 1=y, 2=z), valid after generate().
      //----------------------------------------------------------------------
      double spin(int atom, int component) {
         using namespace quantum::internal;
         // FFT path: read from pre-baked array. generate() guarantees the
         // index is in range, so this is not a silent fallback any more.
         if (sld_fft_n_fine > 0) {
            const size_t idx = static_cast<size_t>(atom) * sld_fft_n_fine + sld_fft_step_index;
            if (component == 0) return sld_fft_spin_x[idx];
            if (component == 1) return sld_fft_spin_y[idx];
            return sld_fft_spin_z[idx];
         }
         // HO path
         if (component == 0) return qn_x_array[atom];
         if (component == 1) return qn_y_array[atom];
         return qn_z_array[atom];
      }

      double phonon(int atom, int component) {
         using namespace quantum::internal;
         if (component == 0) return qn_phonon_x_array[atom];
         if (component == 1) return qn_phonon_y_array[atom];
         return qn_phonon_z_array[atom];
      }

      //----------------------------------------------------------------------
      // Enable per-step export of injected spin noise to a file.
      // Must be called after initialize().
      //----------------------------------------------------------------------
      void enable_spin_noise_export(const std::string& filename, int atom) {
         using namespace quantum::internal;
         sld_export_noise          = true;
         sld_export_noise_filename = filename;
         sld_export_noise_atom     = atom;
      }

      //----------------------------------------------------------------------
      // Append one row to the export file: t  spin_x  spin_y  spin_z
      // No-op when export is not enabled.
      //----------------------------------------------------------------------
      void export_spin_noise_step(double t) {
         using namespace quantum::internal;
         if (!sld_export_noise) return;

         static std::ofstream f;
         if (!f.is_open()) {
            f.open(sld_export_noise_filename, std::ios::trunc);
            if (!f.is_open()) {
               std::cerr << "Warning: cannot open spin-noise export file '"
                         << sld_export_noise_filename << "'\n";
               sld_export_noise = false;
               return;
            }
         }

         const int a = sld_export_noise_atom;
         f << t << '\t'
           << spin(a, 0) << '\t'
           << spin(a, 1) << '\t'
           << spin(a, 2) << '\n';
      }

   } // end of sld_noise namespace

} // end of quantum namespace
