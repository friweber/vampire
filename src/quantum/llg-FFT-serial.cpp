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
//   Serial RK4 time stepper for the FFT noise method.
//
//   Integrates the 9-component state vector y = (S, q, p) per atom using
//   the classical four-stage Runge-Kutta scheme.  Colored noise from the
//   FFT pipeline enters through the effective field H via get_field()
//   (one sample per RK4 step, reused across K1-K4 — LSF_RK4 / HO convention).
//
//   Stage pattern:
//     K1 at t        → y_pred = y_in + (dt/2) * K1
//     K2 at t + dt/2 → y_pred = y_in + (dt/2) * K2
//     K3 at t + dt/2 → y_pred = y_in + dt     * K3
//     K4 at t + dt   → y_final = y_in + (dt/6)(K1 + 2*K2 + 2*K3 + K4)
//
//   Per-atom RK4 primitives (save_initial_state, collect_H_FFT,
//   predict_and_renorm, ...) live in RK4.cpp; this file only orchestrates
//   the RK4 stage sequence.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "quantum.hpp"

// Module headers
#include "internal.hpp"

namespace quantum{
namespace internal{

   //------------------------------------------------------------------------
   // Serial FFT-based LLG step function
   //------------------------------------------------------------------------
   void llg_FFT(){

      // Regenerate the noise window if needed (windowed mode only)
      update_noise_if_needed();

      const int num_atoms = atoms::num_atoms;

      const double dt        = mp::dt;
      const double half_dt   = 0.5 * dt;
      const double dt_over_6 = dt / 6.0;

      double H[3];

      // Snapshot the initial state and evaluate the fields at t.
      for (int atom = 0; atom < num_atoms; ++atom) save_initial_state(atom);

      sim::calculate_spin_fields(0, num_atoms);
      sim::calculate_external_fields(0, num_atoms);

      // Draw the FFT noise once per RK4 step into qn_*_array. The same
      // value is used in every stage's H assembly (LSF_RK4 convention,
      // identical to the HO path).
      draw_noise_all_atoms_FFT(0, num_atoms);

      //=====================================================================
      // K1 Stage (evaluate at t)
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_in_storage[atom].data(), H, k1_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k1_storage[atom].data(), half_dt);
         writeback_spin(atom);
      }
      sim::calculate_spin_fields(0, num_atoms);

      //=====================================================================
      // K2 Stage (evaluate at t + dt/2)
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k2_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k2_storage[atom].data(), half_dt);
         writeback_spin(atom);
      }
      sim::calculate_spin_fields(0, num_atoms);

      //=====================================================================
      // K3 Stage (evaluate at t + dt/2)
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k3_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k3_storage[atom].data(), dt);
         writeback_spin(atom);
      }
      sim::calculate_spin_fields(0, num_atoms);

      //=====================================================================
      // K4 Stage (evaluate at t + dt) + final combine + writeback
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k4_storage[atom].data(), imaterial);
         rk4_combine_and_renorm(atom, dt_over_6);
         writeback_full_state(atom);
      }

      // Advance fine-grained noise time index
      quantum::increment_time();
   }

} // end of internal namespace
} // end of quantum namespace
