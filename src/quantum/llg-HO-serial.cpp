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
//   Serial RK4 time stepper for the HO (harmonic oscillator) noise method.
//
//   Integrates the 9-component state vector y = (S, q, p) per atom.
//   White noise is injected directly into the momentum equation by
//   LL_HO_method() at each RK4 sub-step, so no pre-generated noise
//   array is needed.
//
//   The effective field H here does NOT include noise — it only contains
//   the spin-spin and external fields.
//
//   Per-atom RK4 primitives (save_initial_state, collect_H_HO,
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
   // Serial HO-based LLG step function
   //------------------------------------------------------------------------
   void llg_HO(){

      // Re-fill the T-dependent quantum-noise coefficients if sim::temperature
      // has changed since the previous step (no-op when T is constant).
      refresh_quantum_noise_for_T();

      const int num_atoms = atoms::num_atoms;

      const double dt        = mp::dt;
      const double half_dt   = 0.5 * dt;
      const double dt_over_6 = dt / 6.0;

      double H[3];

      // Snapshot the initial state and evaluate the fields at t.
      for (int atom = 0; atom < num_atoms; ++atom) save_initial_state(atom);

      sim::calculate_spin_fields(0, num_atoms);
      sim::calculate_external_fields(0, num_atoms);

      // Draw thermal noise once per step (held fixed across K1-K4).
      draw_noise_all_atoms_HO(0, num_atoms, dt);

      //=====================================================================
      // K1 Stage (evaluate at t)
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_HO(atom, H);
         LL_HO_method(y_in_storage[atom].data(), H, k1_storage[atom].data(),
                      imaterial, dt,
                      qn_x_array[atom], qn_y_array[atom], qn_z_array[atom]);
         predict_and_renorm(atom, k1_storage[atom].data(), half_dt);
         writeback_spin(atom);
      }
      sim::calculate_spin_fields(0, num_atoms);

      //=====================================================================
      // K2 Stage (evaluate at t + dt/2)
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_HO(atom, H);
         LL_HO_method(y_pred_storage[atom].data(), H, k2_storage[atom].data(),
                      imaterial, dt,
                      qn_x_array[atom], qn_y_array[atom], qn_z_array[atom]);
         predict_and_renorm(atom, k2_storage[atom].data(), half_dt);
         writeback_spin(atom);
      }
      sim::calculate_spin_fields(0, num_atoms);

      //=====================================================================
      // K3 Stage (evaluate at t + dt/2)
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_HO(atom, H);
         LL_HO_method(y_pred_storage[atom].data(), H, k3_storage[atom].data(),
                      imaterial, dt,
                      qn_x_array[atom], qn_y_array[atom], qn_z_array[atom]);
         predict_and_renorm(atom, k3_storage[atom].data(), dt);
         writeback_spin(atom);
      }
      sim::calculate_spin_fields(0, num_atoms);

      //=====================================================================
      // K4 Stage (evaluate at t + dt) + final RK4 combination + writeback
      //=====================================================================
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_HO(atom, H);
         LL_HO_method(y_pred_storage[atom].data(), H, k4_storage[atom].data(),
                      imaterial, dt,
                      qn_x_array[atom], qn_y_array[atom], qn_z_array[atom]);
         rk4_combine_and_renorm(atom, dt_over_6);
         writeback_full_state(atom);
      }

      // Append atom-0 (q_x, q_y, q_z) to the noise file if requested
      if (export_noise) export_ho_noise_step();
   }

} // end of internal namespace
} // end of quantum namespace
