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
//   On-the-fly quantum noise generators for the HO integration method.
//
//   generate_quantum_noise_HO(atom, material, dt)
//
//   Generates quantum-coloured noise for one atom for a single RK4 step.
//   Results are stored in qn_x/y/z_array[atom] and held fixed across all
//   K1-K4 sub-steps.
//
//   Two generators are available, selected by noise_type:
//
//     quantum_zero    -> generate_quantum_noise_orn_uhl()
//       Exact Ornstein-Uhlenbeck integration of n_bath_modes Matsubara modes.
//       Includes a white zero-mode term to reproduce the coth(omega/2T)
//       spectral density (quantum noise with zero-point fluctuations).
//       Unconditionally stable at any time step.
//
//     quantum_no_zero -> generate_quantum_noise_log_bath_opt()
//       Log-spaced auxiliary OU bath with n_bath_modes modes. Lambda range
//       [lambda_lo, lambda_hi] is found at initialisation by a 2-D grid
//       scan that minimises the squared error against the target PSD
//       omega*(coth(omega/2T)-1) (quantum noise without zero-point term).
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// Module headers
#include "internal.hpp"

namespace quantum {
namespace internal {

// Toggle: true = batch-fill noise buffer via std::generate before OU loop; false = generate inline (default)
static bool use_preallocated_noise = true;

   //==========================================================================
   // Exact Ornstein-Uhlenbeck integrator for quantum (zero-point) noise
   //
   // State update (exact): s_n(t+dt) = decay_n * s_n + diffuse_n * xi_n
   //   decay_n   = exp(-nu_n * dt)
   //   diffuse_n = sqrt(1 - decay_n^2)       (Var(s_n) = 1 in steady state)
   //
   // Noise contribution (analytically derived from integrating ds_n over [t,t+dt]):
   //   f_n/dt = sqrt(2T/nu_n)/dt * [(decay_n - 1)*s_n_old + sqrt(1-decay_n^2)*xi_n]
   //          = ou_drift[n]*s_n_old + ou_noise_amp[n]*xi_n
   //
   // White zero-mode term: amp * sqrt(2T/dt) * xi_0  [quantum_zero only]
   //==========================================================================
   void generate_quantum_noise_orn_uhl(const int atom, const double amp, const double /*dt*/) {

      const int    nm   = n_bath_modes;
      const size_t base = static_cast<size_t>(atom) * nm;

      double aux_x = 0.0, aux_y = 0.0, aux_z = 0.0;

      // SoA per-atom buffer: zx[0..nm), zy[nm..2nm), zz[2nm..3nm) — filled once per call
      thread_local static std::vector<double> noise_buf;
      if (use_preallocated_noise) {
         noise_buf.resize(static_cast<size_t>(nm) * 3);
         std::generate(noise_buf.begin(), noise_buf.end(), mtrandom::gaussian);
      }

      // Process two modes per iteration for cache efficiency
      int i = 0;
      for (; i + 1 < nm; i += 2) {

         const double dc0 = ou_decay[i],     dc1 = ou_decay[i + 1];
         const double df0 = ou_diffuse[i],   df1 = ou_diffuse[i + 1];
         const double dr0 = ou_drift[i],     dr1 = ou_drift[i + 1];
         const double na0 = ou_noise_amp[i], na1 = ou_noise_amp[i + 1];

         double& sx0 = matsubara_s_x[base + i];  double& sx1 = matsubara_s_x[base + i + 1];
         double& sy0 = matsubara_s_y[base + i];  double& sy1 = matsubara_s_y[base + i + 1];
         double& sz0 = matsubara_s_z[base + i];  double& sz1 = matsubara_s_z[base + i + 1];

         double zx0, zx1, zy0, zy1, zz0, zz1;
         if (use_preallocated_noise) {
            zx0 = noise_buf[i];              zx1 = noise_buf[i + 1];
            zy0 = noise_buf[nm + i];         zy1 = noise_buf[nm + i + 1];
            zz0 = noise_buf[2 * nm + i];     zz1 = noise_buf[2 * nm + i + 1];
         } else {
            zx0 = mtrandom::gaussian();      zx1 = mtrandom::gaussian();
            zy0 = mtrandom::gaussian();      zy1 = mtrandom::gaussian();
            zz0 = mtrandom::gaussian();      zz1 = mtrandom::gaussian();
         }

         // Noise contribution — uses old state (before update), same xi
         aux_x += dr0 * sx0 + na0 * zx0 + dr1 * sx1 + na1 * zx1;
         aux_y += dr0 * sy0 + na0 * zy0 + dr1 * sy1 + na1 * zy1;
         aux_z += dr0 * sz0 + na0 * zz0 + dr1 * sz1 + na1 * zz1;

         // Exact OU state update (unconditionally stable)
         sx0 = dc0 * sx0 + df0 * zx0;  sx1 = dc1 * sx1 + df1 * zx1;
         sy0 = dc0 * sy0 + df0 * zy0;  sy1 = dc1 * sy1 + df1 * zy1;
         sz0 = dc0 * sz0 + df0 * zz0;  sz1 = dc1 * sz1 + df1 * zz1;
      }

      // Handle odd tail
      if (i < nm) {
         const double dc = ou_decay[i];
         const double df = ou_diffuse[i];
         const double dr = ou_drift[i];
         const double na = ou_noise_amp[i];

         double zx, zy, zz;
         if (use_preallocated_noise) {
            zx = noise_buf[i];  zy = noise_buf[nm + i];  zz = noise_buf[2 * nm + i];
         } else {
            zx = mtrandom::gaussian();  zy = mtrandom::gaussian();  zz = mtrandom::gaussian();
         }

         aux_x += dr * matsubara_s_x[base + i] + na * zx;
         aux_y += dr * matsubara_s_y[base + i] + na * zy;
         aux_z += dr * matsubara_s_z[base + i] + na * zz;

         matsubara_s_x[base + i] = dc * matsubara_s_x[base + i] + df * zx;
         matsubara_s_y[base + i] = dc * matsubara_s_y[base + i] + df * zy;
         matsubara_s_z[base + i] = dc * matsubara_s_z[base + i] + df * zz;
      }

      // White zero-mode: amp * sqrt(2T/dt) * xi
      const double w_x = mtrandom::gaussian();
      const double w_y = mtrandom::gaussian();
      const double w_z = mtrandom::gaussian();

      qn_x_array[atom] = amp * (mats_white_amp * w_x + aux_x);
      qn_y_array[atom] = amp * (mats_white_amp * w_y + aux_y);
      qn_z_array[atom] = amp * (mats_white_amp * w_z + aux_z);
   }

   //==========================================================================
   // Log-spaced OU bath — lambda range optimised by 2-D grid scan
   //
   // K = n_bath_modes auxiliary OU processes with log-spaced rates lambda_k.
   // Range [lambda_lo, lambda_hi] chosen at init by a simple
   // bath_scan_resolution × bath_scan_resolution grid scan over
   // (log_ls, log_le), minimising squared error against the quantum_no_zero
   // target PSD omega*(coth(omega/2T)-1) on a linear ω-grid of
   // bath_scan_omega_points samples up to 15·max(T, ω₀).
   //
   // Algorithm matches the cmp_noise reference (lb_opt_scan).
   //
   // State update (exact OU):
   //   s_k(t+dt) = lb_opt_decay[k] * s_k(t) + lb_opt_coeff_s[k] * xi_k
   //
   // Noise contribution (from updated state):
   //   zeta = amp * sum_k lb_opt_amp[k] * s_k(t+dt)
   //==========================================================================
   void generate_quantum_noise_log_bath_opt(const int atom, const double amp, const double /*dt*/) {

      const int    K    = n_bath_modes;
      const size_t base = static_cast<size_t>(atom) * K;

      double aux_x = 0.0, aux_y = 0.0, aux_z = 0.0;

      // SoA per-atom buffer: zx[0..K), zy[K..2K), zz[2K..3K) — filled once per call
      thread_local static std::vector<double> noise_buf_lb;
      if (use_preallocated_noise) {
         noise_buf_lb.resize(static_cast<size_t>(K) * 3);
         std::generate(noise_buf_lb.begin(), noise_buf_lb.end(), mtrandom::gaussian);
      }

      // Process two modes per iteration for cache efficiency
      int i = 0;
      for (; i + 1 < K; i += 2) {
         const double dc0 = lb_opt_decay[i],   dc1 = lb_opt_decay[i + 1];
         const double cs0 = lb_opt_coeff_s[i], cs1 = lb_opt_coeff_s[i + 1];
         const double la0 = lb_opt_amp[i],     la1 = lb_opt_amp[i + 1];

         double& sx0 = lb_opt_s_x[base + i];     double& sx1 = lb_opt_s_x[base + i + 1];
         double& sy0 = lb_opt_s_y[base + i];     double& sy1 = lb_opt_s_y[base + i + 1];
         double& sz0 = lb_opt_s_z[base + i];     double& sz1 = lb_opt_s_z[base + i + 1];

         double zx0, zx1, zy0, zy1, zz0, zz1;
         if (use_preallocated_noise) {
            zx0 = noise_buf_lb[i];           zx1 = noise_buf_lb[i + 1];
            zy0 = noise_buf_lb[K + i];       zy1 = noise_buf_lb[K + i + 1];
            zz0 = noise_buf_lb[2 * K + i];   zz1 = noise_buf_lb[2 * K + i + 1];
         } else {
            zx0 = mtrandom::gaussian();      zx1 = mtrandom::gaussian();
            zy0 = mtrandom::gaussian();      zy1 = mtrandom::gaussian();
            zz0 = mtrandom::gaussian();      zz1 = mtrandom::gaussian();
         }

         sx0 = dc0 * sx0 + cs0 * zx0;   sx1 = dc1 * sx1 + cs1 * zx1;
         sy0 = dc0 * sy0 + cs0 * zy0;   sy1 = dc1 * sy1 + cs1 * zy1;
         sz0 = dc0 * sz0 + cs0 * zz0;   sz1 = dc1 * sz1 + cs1 * zz1;

         aux_x += la0 * sx0 + la1 * sx1;
         aux_y += la0 * sy0 + la1 * sy1;
         aux_z += la0 * sz0 + la1 * sz1;
      }

      // Handle odd tail
      if (i < K) {
         double zx, zy, zz;
         if (use_preallocated_noise) {
            zx = noise_buf_lb[i];  zy = noise_buf_lb[K + i];  zz = noise_buf_lb[2 * K + i];
         } else {
            zx = mtrandom::gaussian();  zy = mtrandom::gaussian();  zz = mtrandom::gaussian();
         }

         lb_opt_s_x[base + i] = lb_opt_decay[i] * lb_opt_s_x[base + i] + lb_opt_coeff_s[i] * zx;
         lb_opt_s_y[base + i] = lb_opt_decay[i] * lb_opt_s_y[base + i] + lb_opt_coeff_s[i] * zy;
         lb_opt_s_z[base + i] = lb_opt_decay[i] * lb_opt_s_z[base + i] + lb_opt_coeff_s[i] * zz;

         aux_x += lb_opt_amp[i] * lb_opt_s_x[base + i];
         aux_y += lb_opt_amp[i] * lb_opt_s_y[base + i];
         aux_z += lb_opt_amp[i] * lb_opt_s_z[base + i];
      }

      qn_x_array[atom] = amp * aux_x;
      qn_y_array[atom] = amp * aux_y;
      qn_z_array[atom] = amp * aux_z;
   }

   //==========================================================================
   // Cascade OU generator (quantum_no_zero, n_cascade_modes > 1).
   //
   // Each log-bath mode k drives a chain of M=n_cascade_modes OU processes.
   // Level 0 is driven by white noise (standard OU); level j is driven by
   // level j-1. Output is taken from the last level, giving ω^{-2M} roll-off.
   //
   // Discrete update (exact frozen-coefficient Euler, stable for all dt):
   //   s_{k,0}[n+1] = decay_k * s_{k,0}[n] + coeff_s_k * xi
   //   s_{k,j}[n+1] = decay_k * s_{k,j}[n] + casc_coeff_k * s_{k,j-1}[n+1]
   //
   // State layout: lb_opt_casc_s_x[atom * K * M + mode * M + level]
   //==========================================================================
   void generate_quantum_noise_log_bath_opt_cascade(const int atom, const double amp, const double /*dt*/) {

      const int K = n_bath_modes;
      const int M = n_cascade_modes;
      const int base = atom * K * M;

      double aux_x = 0.0, aux_y = 0.0, aux_z = 0.0;

      for (int k = 0; k < K; ++k) {
         const double dk = lb_opt_decay[k];
         const double cs = lb_opt_coeff_s[k];
         const double cc = lb_opt_casc_coeff[k];
         const int idx0 = base + k * M;

         // Level 0: standard OU driven by white noise
         lb_opt_casc_s_x[idx0] = dk * lb_opt_casc_s_x[idx0] + cs * mtrandom::gaussian();
         lb_opt_casc_s_y[idx0] = dk * lb_opt_casc_s_y[idx0] + cs * mtrandom::gaussian();
         lb_opt_casc_s_z[idx0] = dk * lb_opt_casc_s_z[idx0] + cs * mtrandom::gaussian();

         // Levels 1..M-1: each driven by the already-updated level above
         for (int j = 1; j < M; ++j) {
            const int idj = idx0 + j;
            lb_opt_casc_s_x[idj] = dk * lb_opt_casc_s_x[idj] + cc * lb_opt_casc_s_x[idj - 1];
            lb_opt_casc_s_y[idj] = dk * lb_opt_casc_s_y[idj] + cc * lb_opt_casc_s_y[idj - 1];
            lb_opt_casc_s_z[idj] = dk * lb_opt_casc_s_z[idj] + cc * lb_opt_casc_s_z[idj - 1];
         }

         // Accumulate output from last cascade level
         const int idM = idx0 + M - 1;
         aux_x += lb_opt_amp[k] * lb_opt_casc_s_x[idM];
         aux_y += lb_opt_amp[k] * lb_opt_casc_s_y[idM];
         aux_z += lb_opt_amp[k] * lb_opt_casc_s_z[idM];
      }

      qn_x_array[atom] = amp * aux_x;
      qn_y_array[atom] = amp * aux_y;
      qn_z_array[atom] = amp * aux_z;
   }

   //==========================================================================
   // Butterworth low-pass post-filter helpers.
   //
   // Design rationale
   // ─────────────────
   // The log-bath generator produces noise with a 1/ω² spectral envelope
   // at high frequencies, which can drive numerical stiffness.  Rather than
   // modifying the bath itself (the cascade approach), these functions apply
   // a purely post-hoc digital filter to the noise field after it has been
   // written to qn_x/y/z_array.
   //
   // Filter: 4th-order Butterworth low-pass, realised as two cascaded
   // 2nd-order (biquad) IIR sections in Transposed Direct Form II.
   //
   // The two Q factors arise from the Butterworth pole positions:
   //   poles at angles (2k+1)π/(2n) for n=4, k=0..3, in the left half-plane
   //   → two complex-conjugate pairs with damping sin(π/8) and sin(3π/8)
   //   → Q₀ = 1/(2 sin π/8) ≈ 1.3066  (lightly damped, wider pair)
   //      Q₁ = 1/(2 sin 3π/8) ≈ 0.5412 (heavily damped, narrower pair)
   //==========================================================================

   //--------------------------------------------------------------------------
   // calculate_biquad_coeffs
   //
   // Computes the bilinear-transform coefficients for a single 2nd-order
   // Butterworth LP section with quality factor Q and cutoff cutoff_Hz.
   //
   // Pre-warping: ω_a = (2/T) tan(π f_c T)  →  k = tan(π f_c T)
   // This ensures the digital -3 dB point sits exactly at cutoff_Hz.
   //
   // Transfer function after bilinear transform (normalised a0 = 1):
   //
   //         k²          2k²          k²
   //   H = ──────  +  ──────── z⁻¹ + ────── z⁻²
   //         a0           a0            a0
   //       ─────────────────────────────────────────────
   //        1  +  2(k²-1)/a0 z⁻¹  +  (1-k/Q+k²)/a0 z⁻²
   //
   // where  a0 = 1 + k/Q + k²
   //--------------------------------------------------------------------------
   BiquadCoeffs calculate_biquad_coeffs(const double dt_s,
                                         const double cutoff_Hz,
                                         const double Q)
   {
      // Pre-warped analog frequency: tan(π f_c / f_s)
      const double k  = std::tan(M_PI * cutoff_Hz * dt_s);
      const double k2 = k * k;
      const double a0 = 1.0 + k / Q + k2;   // shared normalisation denominator

      BiquadCoeffs c;
      c.b0 =  k2 / a0;                        // numerator: symmetric LP response
      c.b1 =  2.0 * k2 / a0;                  //            b1 = 2*b0
      c.b2 =  k2 / a0;                        //            b2 = b0
      c.a1 =  2.0 * (k2 - 1.0) / a0;          // denominator: resonance position
      c.a2 =  (1.0 - k / Q + k2) / a0;        // denominator: pole damping
      return c;
   }

   //--------------------------------------------------------------------------
   // setup_log_bath_filtered
   //
   // Replaces a direct call to setup_log_bath_opt when the Butterworth filter
   // path is selected.  Delegates bath initialisation to setup_log_bath_opt,
   // then allocates per-atom filter state arrays (all zeroed).
   //
   // Filter coefficients are *not* computed here because dt is not yet
   // available (it appears only as a parameter of the generator).  They are
   // computed lazily on the first call to generate_quantum_noise_log_bath_opt_filtered
   // and cached in butter_stage[].
   //
   // Memory layout of filter_state_x/y/z:
   //   index = atom * 4  +  stage * 2  +  delay
   //   stage  ∈ {0, 1}  — biquad section index
   //   delay  ∈ {0, 1}  — s1, s2 of the Transposed DFI signal flow graph
   //--------------------------------------------------------------------------
   void setup_log_bath_filtered(const int num_atoms_total,
                                  const double cutoff_frequency_THz)
   {
      // Initialise the underlying log-bath (λ-scan, OU mode allocation, etc.)
      setup_log_bath_opt(num_atoms_total);

      // Store cutoff; coefficients are computed at first generate call (needs dt)
      butter_cutoff_Hz    = cutoff_frequency_THz * 1.0e12;
      butter_coeffs_valid = false;

      // Allocate flat per-atom state: 2 sections × 2 delays = 4 doubles per atom
      const std::size_t n_states = static_cast<std::size_t>(num_atoms_total) * 4;
      filter_state_x.assign(n_states, 0.0);
      filter_state_y.assign(n_states, 0.0);
      filter_state_z.assign(n_states, 0.0);
   }

   //--------------------------------------------------------------------------
   // generate_quantum_noise_log_bath_opt_filtered
   //
   // Wrapper that:
   //   1. Calls generate_quantum_noise_log_bath_opt (writes raw quantum noise
   //      into qn_x/y/z_array[atom]).
   //   2. Passes the raw field through two cascaded Butterworth biquad stages
   //      using per-atom Transposed Direct Form II state buffers.
   //   3. Overwrites qn_x/y/z_array[atom] with the filtered field.
   //
   // Transposed Direct Form II update for one biquad section:
   //
   //   y[n]   =  b0·x[n] + s1            ← output; uses OLD s1
   //   s1[n]  =  b1·x[n] − a1·y[n] + s2  ← new s1; uses OLD s2
   //   s2[n]  =  b2·x[n] − a2·y[n]       ← new s2
   //
   // This form is numerically superior to Direct Form I (no subtraction of
   // large nearly-equal numbers) and requires only 2 delay elements per section.
   //--------------------------------------------------------------------------
   void generate_quantum_noise_log_bath_opt_filtered(const int atom,
                                                      const double amp,
                                                      const double dt)
   {
      // ── Lazy coefficient initialisation (runs exactly once per run) ──────────
      // dt here is mp::dt (VAMPIRE normalised: dt_SI × γ_SI).
      // Convert to physical seconds: dt_s = dt / γ_SI.
      if (!butter_coeffs_valid) {
         constexpr double gamma_SI = 1.76085963e11;  // rad/(s·T)
         const double dt_s = dt / gamma_SI;

         // 4th-order Butterworth pole Q factors (exact analytic values):
         //   Q₀ = 1/(2 sin(π/8))  — from the pole pair at angle 3π/8
         //   Q₁ = 1/(2 sin(3π/8)) — from the pole pair at angle 5π/8
         constexpr double Q0 = 1.0 / (2.0 * 0.38268343236508977);  // ≈ 1.3066
         constexpr double Q1 = 1.0 / (2.0 * 0.92387953251128676);  // ≈ 0.5412
         butter_stage[0] = calculate_biquad_coeffs(dt_s, butter_cutoff_Hz, Q0);
         butter_stage[1] = calculate_biquad_coeffs(dt_s, butter_cutoff_Hz, Q1);
         butter_coeffs_valid = true;
      }

      // ── Step 1: generate raw quantum-coloured noise ───────────────────────────
      generate_quantum_noise_log_bath_opt(atom, amp, dt);

      // ── Step 2: apply cascaded biquad filter to each Cartesian axis ───────────
      // The two stages are applied in series; the output of stage 0 is the
      // input to stage 1 (in-place in the local x_x / x_y / x_z variables).
      double x_x = qn_x_array[atom];
      double x_y = qn_y_array[atom];
      double x_z = qn_z_array[atom];

      const int base = atom * 4;   // first state index for this atom

      for (int stage = 0; stage < 2; ++stage) {
         const int         idx = base + stage * 2;  // s1 at idx, s2 at idx+1
         const BiquadCoeffs& c = butter_stage[stage];

         // — x component ───────────────────────────────────────────────────────
         {
            // Read old state (must happen before any writes)
            const double s1_old = filter_state_x[idx];
            const double s2_old = filter_state_x[idx + 1];
            const double y  = c.b0 * x_x + s1_old;             // output
            filter_state_x[idx]     = c.b1 * x_x - c.a1 * y + s2_old;
            filter_state_x[idx + 1] = c.b2 * x_x - c.a2 * y;
            x_x = y;
         }

         // — y component ───────────────────────────────────────────────────────
         {
            const double s1_old = filter_state_y[idx];
            const double s2_old = filter_state_y[idx + 1];
            const double y  = c.b0 * x_y + s1_old;
            filter_state_y[idx]     = c.b1 * x_y - c.a1 * y + s2_old;
            filter_state_y[idx + 1] = c.b2 * x_y - c.a2 * y;
            x_y = y;
         }

         // — z component ───────────────────────────────────────────────────────
         {
            const double s1_old = filter_state_z[idx];
            const double s2_old = filter_state_z[idx + 1];
            const double y  = c.b0 * x_z + s1_old;
            filter_state_z[idx]     = c.b1 * x_z - c.a1 * y + s2_old;
            filter_state_z[idx + 1] = c.b2 * x_z - c.a2 * y;
            x_z = y;
         }
      }

      // ── Step 3: overwrite global noise arrays with filtered result ────────────
      qn_x_array[atom] = x_x;
      qn_y_array[atom] = x_y;
      qn_z_array[atom] = x_z;
   }

   //==========================================================================
   // Dispatcher — routes to the integrator selected by noise_type
   //==========================================================================
   void generate_quantum_noise_HO(const int atom, const int material, const double dt) {
      const double amp = std::sqrt(  material_gamma_array[material]
                                   * material_A_array[material]
                                   / material_S0_array[material]);
      switch (noise_type) {
         case quantum_zero:    generate_quantum_noise_orn_uhl(atom, amp, dt);        break;
         case quantum_no_zero:
            if (butter_cutoff_Hz > 0.0)
               generate_quantum_noise_log_bath_opt_filtered(atom, amp, dt);
            else if (n_cascade_modes > 1)
               generate_quantum_noise_log_bath_opt_cascade(atom, amp, dt);
            else
               generate_quantum_noise_log_bath_opt(atom, amp, dt);
            break;
         default: break;
      }
   }

   //==========================================================================
   // Refresh: write the T-dependent OU coefficient arrays at a given T.
   //
   //   ν_i      = 2π(i+1) T_scaled
   //   decay_i  = exp(-ν_i dt)
   //   diffuse_i= sqrt(1 - decay_i²)
   //   drift_i  = sqrt(2T/ν_i) · (decay_i − 1) / dt
   //   amp_i    = sqrt(2T/ν_i) · diffuse_i / dt
   //   white_amp= sqrt(2T/dt)
   //
   // Called on every T change.  Cheap: O(n_bath_modes).
   //==========================================================================
   void refresh_orn_uhl_coeffs(const double T_scaled, const double dt) {
      const double inv_dt = 1.0 / dt;
      for (int i = 0; i < n_bath_modes; ++i) {
         const double nu        = 2.0 * M_PI * (i + 1) * T_scaled;
         const double decay_i   = std::exp(-nu * dt);
         const double diffuse_i = std::sqrt(1.0 - decay_i * decay_i);
         ou_decay[i]     = decay_i;
         ou_diffuse[i]   = diffuse_i;
         ou_drift[i]     = std::sqrt(2.0 * T_scaled / nu) * (decay_i - 1.0) * inv_dt;
         ou_noise_amp[i] = std::sqrt(2.0 * T_scaled / nu) * diffuse_i * inv_dt;
      }
      mats_white_amp = std::sqrt(2.0 * T_scaled * inv_dt);
   }

   //==========================================================================
   // Setup: exact Ornstein-Uhlenbeck path (quantum_zero noise type)
   //
   // One-time allocation of per-atom Matsubara state and coefficient arrays.
   // The T-dependent coefficient values are filled by refresh_orn_uhl_coeffs.
   //==========================================================================
   void setup_orn_uhl(const int num_atoms_total) {

      const size_t matsubara_size = static_cast<size_t>(num_atoms_total) * n_bath_modes;
      matsubara_s_x.assign(matsubara_size, 0.0);
      matsubara_s_y.assign(matsubara_size, 0.0);
      matsubara_s_z.assign(matsubara_size, 0.0);

      ou_decay.resize(n_bath_modes);
      ou_diffuse.resize(n_bath_modes);
      ou_drift.resize(n_bath_modes);
      ou_noise_amp.resize(n_bath_modes);

      // Initial coefficient fill at the current temperature; later refreshes
      // happen via refresh_quantum_noise_for_T() if sim::temperature changes.
      const double T_scaled = scale_temperature(sim::temperature);
      refresh_orn_uhl_coeffs(T_scaled, mp::dt);
      last_T_scaled = T_scaled;   // cache is now in sync with the coefficients

      std::cout << "  Generator         : orn-uhl  (exact OU, "
                << n_bath_modes << " Matsubara modes)" << std::endl;
   }

   //==========================================================================
   // Objective for the log-bath λ-range search.
   //
   // Returns the total squared error between the model bath PSD
   //    S_F(ω) = Σ_k 2 λ_k² G_k dl_k / (ω² + λ_k²)
   // and the target
   //    S_F_tgt(ω) = ω · (coth(ω/2T) − 1)
   // for n_bath_modes modes log-spaced in [10^log_ls, 10^log_le].
   // dl_k uses forward differences (k=0..K-2 sum).
   //==========================================================================
   double lb_objective(const double log_ls, const double log_le,
                       const double T_scaled, const int K,
                       const std::vector<double>& scan_om) {
      if (log_le <= log_ls + 0.05) {
         return std::numeric_limits<double>::infinity();
      }
      std::vector<double> lam(K);
      for (int k = 0; k < K; ++k)
         lam[k] = std::pow(10.0, log_ls + static_cast<double>(k) / (K - 1) * (log_le - log_ls));

      double tot_err = 0.0;
      const int N_OM = static_cast<int>(scan_om.size());
      for (int ii = 0; ii < N_OM; ++ii) {
         const double om = scan_om[ii];
         double S_model = 0.0;
         for (int k = 0; k < K - 1; ++k) {
            const double dl = lam[k + 1] - lam[k];
            const double x  = lam[k] / (2.0 * T_scaled);
            const double ct = (std::fabs(x) < 1e-10) ? 1.0 / x : 1.0 / std::tanh(x);
            const double G  = ct - 1.0;
            S_model += 2.0 * lam[k] * lam[k] * G * dl / (om * om + lam[k] * lam[k]);
         }
         const double x_om  = om / (2.0 * T_scaled);
         const double ct_om = (std::fabs(x_om) < 1e-10) ? 1.0 / x_om : 1.0 / std::tanh(x_om);
         const double S_tgt = om * (ct_om - 1.0);
         const double diff  = S_model - S_tgt;
         tot_err += diff * diff;
      }
      return tot_err;
   }

   //==========================================================================
   // Refresh: write the T-dependent log-bath coefficient arrays at a given T.
   //
   // The λ-range (lb_opt_lambdas) is FIXED at init time — re-running the
   // grid scan per T change would dominate the per-step cost.  Only the
   // per-mode amplitudes and OU decay/diffuse coefficients are refreshed:
   //
   //   raw_amp_k = sqrt(λ_k · (coth(λ_k/2T) − 1) · dλ_k)
   //   renorm    = sqrt( target(ω_ref) / Σ_k raw_amp_k² · 2λ_k / (ω_ref²+λ_k²) )
   //   amp_k     = raw_amp_k · renorm
   //   decay_k   = exp(−λ_k · dt)
   //   coeff_s_k = sqrt(1 − decay_k²)
   //
   // Note: decay_k and coeff_s_k depend on dt only, not on T.  Recomputing
   // them is still cheap (a few `exp` / `sqrt` per mode) and keeps the
   // refresh self-contained.
   //==========================================================================
   void refresh_log_bath_opt_coeffs(const double T_scaled, const double dt) {

      const int K = n_bath_modes;
      const double* lam = lb_opt_lambdas.data();

      // Pass 1: raw amplitudes
      std::vector<double> raw_amp(K, 0.0);
      for (int k = 0; k < K; ++k) {
         const double dl = (k < K - 1) ? lam[k + 1] - lam[k]
                                       : lam[k]     - lam[k - 1];
         const double x  = lam[k] / (2.0 * T_scaled);
         const double ct = (std::fabs(x) < 1e-10) ? 1.0 / x : 1.0 / std::tanh(x);
         const double G  = ct - 1.0;
         raw_amp[k] = (G > 0.0) ? std::sqrt(lam[k] * G * dl) : 0.0;
      }
      // Pass 2: DC renorm at ω_ref << T_scaled
      const double omega_ref = T_scaled * 1e-3;
      double sum_L = 0.0;
      for (int k = 0; k < K; ++k)
         sum_L += raw_amp[k] * raw_amp[k]
                  * 2.0 * lam[k] / (omega_ref * omega_ref + lam[k] * lam[k]);
      const double x_ref  = omega_ref / (2.0 * T_scaled);
      const double target = omega_ref * (1.0 / std::tanh(x_ref) - 1.0);
      const double renorm = (sum_L > 0.0) ? std::sqrt(target / sum_L) : 1.0;
      // Pass 3: write final arrays
      for (int k = 0; k < K; ++k) {
         const double d = std::exp(-lam[k] * dt);
         lb_opt_decay[k]   = d;
         lb_opt_coeff_s[k] = std::sqrt(1.0 - d * d);
         lb_opt_amp[k]     = raw_amp[k] * renorm;
         if (!lb_opt_casc_coeff.empty())
            // Coupling for the cascade chain.  Using λ^{1/M} (M-th root) ensures
            // that the low-frequency output power degrades only as λ^{-2(M-1)/M}
            // per mode rather than λ^{-(M-1)} with the fixed sqrt, so the cascade
            // output does not collapse to zero for deeper chains.
            // For M=2: λ^{1/2} = sqrt(λ) — identical to the original formula.
            lb_opt_casc_coeff[k] = (1.0 - d)
                                 / std::pow(lam[k], 1.0 / n_cascade_modes);
      }
   }

   //==========================================================================
   // Setup: opt-range log-spaced OU bath (quantum_no_zero noise type)
   //
   // One-time work:
   //   * allocate per-atom state arrays
   //   * run the 2-D (log_ls, log_le) grid scan to choose the λ-range
   //     (algorithm and defaults match cmp_noise.cpp's lb_opt_scan)
   //   * cache λ-range in lb_opt_lambdas (fixed across T changes)
   //
   // The T-dependent coefficient arrays (decay / coeff_s / amp) are then
   // filled by refresh_log_bath_opt_coeffs at the current temperature.
   //==========================================================================
   void setup_log_bath_opt(const int num_atoms_total) {

      const double T_scaled = scale_temperature(sim::temperature);
      const double dt       = mp::dt;

      const size_t lb_opt_size = static_cast<size_t>(num_atoms_total) * n_bath_modes;
      lb_opt_s_x.assign(lb_opt_size, 0.0);
      lb_opt_s_y.assign(lb_opt_size, 0.0);
      lb_opt_s_z.assign(lb_opt_size, 0.0);

      const size_t casc_size = static_cast<size_t>(num_atoms_total) * n_bath_modes * n_cascade_modes;
      lb_opt_casc_s_x.assign(casc_size, 0.0);
      lb_opt_casc_s_y.assign(casc_size, 0.0);
      lb_opt_casc_s_z.assign(casc_size, 0.0);
      lb_opt_casc_coeff.resize(n_bath_modes);

      const double omega0_opt = material_omega0_array.empty() ? T_scaled : material_omega0_array[0];

      // Frequency grid the objective is summed over.
      //
      // Linear from ω_max/N to ω_max with ω_max = 15·max(T_scaled, ω₀),
      // matching the cmp_noise reference algorithm.  The linear grid
      // concentrates samples around the Lorentzian peak where the spin
      // actually couples; extending only to 15× max keeps the optimiser
      // from being dragged by the unavoidable 1/ω² polynomial tail of the
      // finite Lorentzian sum at frequencies the physics never sees.
      //
      // Point count is user-tunable via quantum:bath-scan-omega-points
      // (default 500 = cmp_noise reference).
      const int    N_OM   = bath_scan_omega_points;
      const double om_max = std::max(T_scaled, omega0_opt) * 15.0;
      std::vector<double> scan_om(N_OM);
      for (int ii = 0; ii < N_OM; ++ii)
         scan_om[ii] = (ii + 1) * om_max / N_OM;

      // Bounding box for the 2-D grid scan.  The asymmetric inner offsets
      // (+1.0 and -0.5) are physics-motivated and constant: λ_lo should sit
      // below T, λ_hi above max(T, ω₀).  The outer offsets are widened/
      // narrowed by the user-tunable quantum:bath-scan-decades (default 3.5).
      const double D      = bath_scan_decades;
      const double ls_min = std::log10(T_scaled)                       - D;
      const double ls_max = std::log10(T_scaled)                       + 1.0;
      const double le_min = std::log10(std::max(T_scaled, omega0_opt)) - 0.5;
      const double le_max = std::log10(std::max(T_scaled, omega0_opt)) + D;

      // Simple 2-D grid scan over (log_ls, log_le) — matches cmp_noise's
      // lb_opt_scan with resolution configurable via
      // quantum:bath-scan-resolution (default 50).  At the default the
      // total cost is 50² × N_OM ≈ 1.25 M cheap evals — sub-second.
      const int n_scan = bath_scan_resolution;
      double best_err = std::numeric_limits<double>::infinity();
      double best_ls  = ls_min;
      double best_le  = le_min;
      for (int i = 0; i < n_scan; ++i) {
         const double log_ls = ls_min + i * (ls_max - ls_min) / (n_scan - 1);
         for (int j = 0; j < n_scan; ++j) {
            const double log_le = le_min + j * (le_max - le_min) / (n_scan - 1);
            if (log_le <= log_ls + 0.05) continue;
            const double err = lb_objective(log_ls, log_le, T_scaled, n_bath_modes, scan_om);
            if (err < best_err) {
               best_err = err;
               best_ls  = log_ls;
               best_le  = log_le;
            }
         }
      }

      // Cache the λ-range for the per-T refresh path
      lb_opt_lambdas.resize(n_bath_modes);
      for (int k = 0; k < n_bath_modes; ++k)
         lb_opt_lambdas[k] = std::pow(10.0, best_ls
                                  + static_cast<double>(k) / (n_bath_modes - 1)
                                  * (best_le - best_ls));

      // Initial coefficient fill at the current temperature
      lb_opt_decay.resize(n_bath_modes);
      lb_opt_coeff_s.resize(n_bath_modes);
      lb_opt_amp.resize(n_bath_modes);
      refresh_log_bath_opt_coeffs(T_scaled, dt);
      last_T_scaled = T_scaled;   // cache is now in sync with the coefficients

      std::cout << "  Generator         : log-bath-opt-range  ("
                << n_bath_modes << " log-spaced OU modes)" << std::endl;
      std::cout << "                       lambda = [" << std::pow(10.0, best_ls)
                << ", " << std::pow(10.0, best_le) << "]  ("
                << n_scan << "x" << n_scan << " grid scan, err=" << best_err << ")"
                << std::endl;
   }

   //==========================================================================
   // Top-level dynamic-T refresh hook.
   //
   // Called at the top of every LLG step (HO method).  Compares the current
   // scale_temperature(sim::temperature) against the cached last_T_scaled:
   //   * if equal (constant-T runs), this is one float compare and returns.
   //   * if different (temperature pulse, laser pulse, …), re-runs the
   //     appropriate refresh helper and updates last_T_scaled.
   //==========================================================================
   void refresh_quantum_noise_for_T() {
      const double T_now = scale_temperature(sim::temperature);
      if (T_now == last_T_scaled) return;
      const double dt = mp::dt;
      if (noise_type == quantum_zero) {
         refresh_orn_uhl_coeffs(T_now, dt);
      } else if (noise_type == quantum_no_zero) {
         refresh_log_bath_opt_coeffs(T_now, dt);
      }
      last_T_scaled = T_now;
   }

   //==========================================================================
   // HO noise export — append one row per LLG step (rank 0 only)
   //
   // For the HO method, the "noise" that drives the spin is the auxiliary
   // oscillator coordinate q (dS/dt = S × (H + q)). Saving q(t) for atom 0
   // is the direct analogue of what export_noise_data() writes for the FFT
   // path. Called from llg_HO / llg_HO_mpi at end-of-step, when export_noise
   // is true.
   //==========================================================================
   void export_ho_noise_step() {
      #ifdef MPICF
      if (vmpi::my_rank != 0) return;
      #endif

      // Skip the first noise_export_burnin_steps calls so the recorded
      // trace starts in the HO's statistical steady state.  Mirrors
      // cmp_noise.cpp's n_burnin = min(nsteps/5, 20000) — without it,
      // the q-from-zero ramp dominates the low-frequency PSD and the
      // shape no longer matches the Lorentzian target.
      const uint64_t step_now = noise_export_step_count++;
      if (step_now < noise_export_burnin_steps) return;

      // First post-burn-in call: truncate + write header.
      // Subsequent calls: append.
      const std::ios_base::openmode mode = noise_export_header_written
         ? (std::ios_base::out | std::ios_base::app)
         : std::ios_base::out;

      std::ofstream out(export_noise_filename, mode);
      if (!out.is_open()) {
         std::cerr << "Error: cannot open HO noise export file "
                   << export_noise_filename << std::endl;
         err::vexit();
      }

      if (!noise_export_header_written) {
         out << "# HO quantum noise export — auxiliary oscillator q for atom 0\n";
         out << "# noise_type: " << noise_type_name(noise_type) << "\n";
         out << "# method:     " << llg_method_short_name(llg_method) << "\n";
         out << "# dt:         " << mp::dt << " s\n";
         out << "# Column 1: time (s)\n";
         out << "# Column 2: q_x (atom 0)\n";
         out << "# Column 3: q_y (atom 0)\n";
         out << "# Column 4: q_z (atom 0)\n";
         noise_export_header_written = true;
      }

      const double t = static_cast<double>(sim::time) * mp::dt;
      out << t << " "
          << q_x_array[0] << " "
          << q_y_array[0] << " "
          << q_z_array[0] << "\n";
   }

} // end of internal namespace
} // end of quantum namespace
