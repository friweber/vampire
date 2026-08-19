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
//   Shared RK4 plumbing for the quantum thermostat time-steppers.
//
//   Equations of motion (per atom, 9-component state y = (S, q, p)):
//     LL_HO_method  — HO method: white noise injected into the momentum equation
//     LL_FFT_method — FFT method: noise pre-baked into the effective field H
//
//   Per-atom RK4 primitives used by all four time-steppers:
//     save_initial_state     — atoms arrays → y_in_storage
//     collect_H_HO           — H = spin_field + external_field
//     collect_H_FFT          — H = spin_field + external_field + qn_{x,y,z}_array
//     predict_and_renorm     — y_pred = y_in + coeff*k; renormalise spin
//     rk4_combine_and_renorm — y_pred = y_in + dt/6*(k1+2k2+2k3+k4); renorm spin
//     writeback_spin         — y_pred[0..2] → atoms spin arrays
//     writeback_full_state   — y_pred[0..8] → atoms spin + q + p arrays
//     renormalize_spin       — utility used by the predict/combine helpers
//
//   Noise draw (HO method only — FFT noise is pre-generated):
//     draw_noise_all_atoms_HO  — fills qn_{x,y,z}_array[lo..hi) (HO)
//     draw_noise_all_atoms_FFT — fills qn_{x,y,z}_array[lo..hi) (FFT)
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"

// Module headers
#include "internal.hpp"

namespace quantum{
   namespace internal{

      //=====================================================================
      // HO method: equations of motion with on-the-fly white noise
      //
      // The noise prefactor is: sqrt(2 * Gamma * A * T / (S0 * dt))
      // ensuring correct fluctuation-dissipation balance.
      //=====================================================================
      void LL_HO_method(const double* y, const double* H, double* dydt,
                        const int material, const double dt,
                        const double noise_x, const double noise_y, const double noise_z) {

         const double A      = material_A_array[material];
         const double Gamma  = material_gamma_array[material];
         const double omega0 = material_omega0_array[material];

         // Noise is pre-drawn once per full RK4 step and passed in,
         // ensuring the same noise realization is used across K1-K4.

         // dS/dt = S x (H + q)
         dydt[0] = y[1]*(H[2] + y[5]) - y[2]*(H[1] + y[4]);
         dydt[1] = y[2]*(H[0] + y[3]) - y[0]*(H[2] + y[5]);
         dydt[2] = y[0]*(H[1] + y[4]) - y[1]*(H[0] + y[3]);

         // dq/dt = p
         dydt[3] = y[6];
         dydt[4] = y[7];
         dydt[5] = y[8];

         // dp/dt = -omega0^2 * q - Gamma * p + A * S + noise
         dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0] + noise_x;
         dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1] + noise_y;
         dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2] + noise_z;
      }

      //=====================================================================
      // FFT method: equations of motion (noise is in H, not momentum)
      //
      // The effective field H passed to this function already includes the
      // pre-generated colored noise from the FFT pipeline.
      //=====================================================================
      void LL_FFT_method(const double* y, const double* H, double* dydt,
                         const int material) {

         const double A      = material_A_array[material];
         const double Gamma  = material_gamma_array[material];
         const double omega0 = material_omega0_array[material];

         // dS/dt = S x (H + q)
         dydt[0] = y[1]*(H[2] + y[5]) - y[2]*(H[1] + y[4]);
         dydt[1] = y[2]*(H[0] + y[3]) - y[0]*(H[2] + y[5]);
         dydt[2] = y[0]*(H[1] + y[4]) - y[1]*(H[0] + y[3]);

         // dq/dt = p
         dydt[3] = y[6];
         dydt[4] = y[7];
         dydt[5] = y[8];

         // dp/dt = -omega0^2 * q - Gamma * p + A * S   (no noise term)
         dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0];
         dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1];
         dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2];
      }

      //=====================================================================
      // Per-atom RK4 primitives
      //=====================================================================

      // Save current (S, q, p) state into y_in_storage[atom].
      void save_initial_state(const int atom) {
         y_in_storage[atom][0] = atoms::x_spin_array[atom];
         y_in_storage[atom][1] = atoms::y_spin_array[atom];
         y_in_storage[atom][2] = atoms::z_spin_array[atom];
         y_in_storage[atom][3] = q_x_array[atom];
         y_in_storage[atom][4] = q_y_array[atom];
         y_in_storage[atom][5] = q_z_array[atom];
         y_in_storage[atom][6] = p_x_array[atom];
         y_in_storage[atom][7] = p_y_array[atom];
         y_in_storage[atom][8] = p_z_array[atom];
      }

      // Normalise the spin part (first three components) of a 9-vector in place.
      void renormalize_spin(double* y) {
         const double S_mag   = std::sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
         const double inv_mag = 1.0 / S_mag;
         y[0] *= inv_mag;
         y[1] *= inv_mag;
         y[2] *= inv_mag;
      }

      // Effective field H for the HO method (no noise — it enters via LL_HO_method).
      void collect_H_HO(const int atom, double H[3]) {
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom];
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom];
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom];
      }

      // Effective field H for the FFT method, including pre-generated colored noise.
      // The noise component reads from qn_{x,y,z}_array, populated once per
      // RK4 step by draw_noise_all_atoms_FFT() — same convention as the HO
      // path and as LSF_RK4 (one sample per step, reused across K1-K4).
      void collect_H_FFT(const int atom, double H[3]) {
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom]
              + qn_x_array[atom];
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom]
              + qn_y_array[atom];
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
              + qn_z_array[atom];
      }

      // Sample the FFT noise field for atoms [lo, hi) into qn_{x,y,z}_array.
      // One sample per atom per RK4 step — held fixed across K1-K4.
      void draw_noise_all_atoms_FFT(const int lo, const int hi) {
         for (int atom = lo; atom < hi; ++atom) {
            qn_x_array[atom] = quantum::get_field(atom, 0);
            qn_y_array[atom] = quantum::get_field(atom, 1);
            qn_z_array[atom] = quantum::get_field(atom, 2);
         }
      }

      // y_pred_storage[atom] = y_in_storage[atom] + coeff * k_stage,
      // then renormalise the spin part.
      void predict_and_renorm(const int atom, const double* k_stage, const double coeff) {
         double*       yp = y_pred_storage[atom].data();
         const double* yi = y_in_storage[atom].data();
         for (size_t i = 0; i < 9; ++i) {
            yp[i] = yi[i] + coeff * k_stage[i];
         }
         renormalize_spin(yp);
      }

      // y_pred = y_in + (dt/6)(k1 + 2*k2 + 2*k3 + k4), then renormalise spin.
      // Does NOT write back to the atoms arrays — caller decides when.
      void rk4_combine_and_renorm(const int atom, const double dt_over_6) {
         double*       yp = y_pred_storage[atom].data();
         const double* yi = y_in_storage[atom].data();
         const double* k1 = k1_storage[atom].data();
         const double* k2 = k2_storage[atom].data();
         const double* k3 = k3_storage[atom].data();
         const double* k4 = k4_storage[atom].data();
         for (size_t i = 0; i < 9; ++i) {
            yp[i] = yi[i] + dt_over_6 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
         }
         renormalize_spin(yp);
      }

      // Copy spin components of y_pred back to atoms::*_spin_array.
      // Used between RK4 stages so that spin-field recalculation sees the
      // current predictor.
      void writeback_spin(const int atom) {
         atoms::x_spin_array[atom] = y_pred_storage[atom][0];
         atoms::y_spin_array[atom] = y_pred_storage[atom][1];
         atoms::z_spin_array[atom] = y_pred_storage[atom][2];
      }

      // Copy all 9 components of y_pred back to the per-atom state arrays.
      // Used after the final RK4 combination to commit the new state.
      void writeback_full_state(const int atom) {
         atoms::x_spin_array[atom] = y_pred_storage[atom][0];
         atoms::y_spin_array[atom] = y_pred_storage[atom][1];
         atoms::z_spin_array[atom] = y_pred_storage[atom][2];
         q_x_array[atom] = y_pred_storage[atom][3];
         q_y_array[atom] = y_pred_storage[atom][4];
         q_z_array[atom] = y_pred_storage[atom][5];
         p_x_array[atom] = y_pred_storage[atom][6];
         p_y_array[atom] = y_pred_storage[atom][7];
         p_z_array[atom] = y_pred_storage[atom][8];
      }

      //=====================================================================
      // HO noise draw (one realisation per atom per RK4 step)
      //
      // Classical: white Gaussian scaled by the FDT prefactor
      //              sqrt(2 * Gamma * A * T / (S0 * dt)).
      // Quantum:   delegated to generate_quantum_noise_HO(), which dispatches
      //              on noise_type to orn_uhl (quantum_zero) or
      //              log_bath_opt (quantum_no_zero).
      // Result stored in qn_{x,y,z}_array[atom] and held fixed across K1-K4.
      //=====================================================================
      void draw_noise_all_atoms_HO(const int lo, const int hi, const double dt) {
         const double T_scaled = scale_temperature(sim::temperature);

         for (int atom = lo; atom < hi; ++atom) {
            const int imaterial = atoms::type_array[atom];
            if (noise_type == classical) {
               const double A     = material_A_array[imaterial];
               const double Gamma = material_gamma_array[imaterial];
               const double S0    = material_S0_array[imaterial];
               const double noise_prefactor = std::sqrt(2.0 * Gamma * A * T_scaled / (S0 * dt));
               qn_x_array[atom] = noise_prefactor * mtrandom::gaussian();
               qn_y_array[atom] = noise_prefactor * mtrandom::gaussian();
               qn_z_array[atom] = noise_prefactor * mtrandom::gaussian();
            } else {
               generate_quantum_noise_HO(atom, imaterial, dt);
            }
         }
      }

   } // end of internal namespace
} // end of quantum namespace
