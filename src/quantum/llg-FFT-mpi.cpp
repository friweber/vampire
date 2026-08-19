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
//   MPI-parallel RK4 time stepper for the FFT noise method.
//
//   Same RK4 scheme as the serial version but with halo-swap communication
//   overlapped with core-atom computation at each stage:
//
//     1. Initiate non-blocking halo swap
//     2. Compute K_i for core atoms (no dependency on boundary)
//     3. Complete halo swap
//     4. Compute K_i for boundary atoms
//     5. Update spin arrays and proceed to next stage
//
//   Colored noise from the FFT pipeline is drawn once per RK4 step into
//   qn_{x,y,z}_array (matching the HO / LSF_RK4 convention) and reused
//   across all four sub-stages via collect_H_FFT().
//
//   Per-atom RK4 primitives (save_initial_state, collect_H_FFT,
//   predict_and_renorm, ...) live in RK4.cpp; this file only orchestrates
//   the RK4 stage sequence + MPI synchronisation.
//
//------------------------------------------------------------------------------

#ifdef MPICF

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "quantum.hpp"
#include "vmpi.hpp"

// Module headers
#include "internal.hpp"

namespace quantum{
namespace internal{

   //------------------------------------------------------------------------
   // MPI FFT-based LLG step function
   //------------------------------------------------------------------------
   void llg_FFT_mpi(){

      // Regenerate the noise window if needed (windowed mode only)
      update_noise_if_needed();

      // MPI atom ranges
      const int pre_comm_si  = 0;
      const int pre_comm_ei  = vmpi::num_core_atoms;
      const int post_comm_si = vmpi::num_core_atoms;
      const int post_comm_ei = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

      const double dt        = mp::dt;
      const double half_dt   = 0.5 * dt;
      const double dt_over_6 = dt / 6.0;

      double H[3];

      // Initial halo swap to get valid boundary atom data
      vmpi::mpi_init_halo_swap();
      vmpi::mpi_complete_halo_swap();

      // Snapshot the initial state for all atoms (core + boundary)
      for (int atom = pre_comm_si; atom < post_comm_ei; ++atom) save_initial_state(atom);

      // Draw the FFT noise once per RK4 step into qn_*_array, reused across
      // K1-K4 (LSF_RK4 / HO convention).
      draw_noise_all_atoms_FFT(pre_comm_si, post_comm_ei);

      //=====================================================================
      // K1 Stage (evaluate at t)
      //=====================================================================
      vmpi::mpi_init_halo_swap();
      sim::calculate_spin_fields(pre_comm_si, pre_comm_ei);
      sim::calculate_external_fields(pre_comm_si, post_comm_ei);

      for (int atom = pre_comm_si; atom < pre_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_in_storage[atom].data(), H, k1_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k1_storage[atom].data(), half_dt);
      }

      vmpi::mpi_complete_halo_swap();
      sim::calculate_spin_fields(post_comm_si, post_comm_ei);
      sim::calculate_external_fields(post_comm_si, post_comm_ei);

      for (int atom = post_comm_si; atom < post_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_in_storage[atom].data(), H, k1_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k1_storage[atom].data(), half_dt);
      }

      for (int atom = pre_comm_si; atom < post_comm_ei; ++atom) writeback_spin(atom);

      //=====================================================================
      // K2 Stage (evaluate at t + dt/2)
      //=====================================================================
      vmpi::mpi_init_halo_swap();
      sim::calculate_spin_fields(pre_comm_si, pre_comm_ei);

      for (int atom = pre_comm_si; atom < pre_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k2_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k2_storage[atom].data(), half_dt);
      }

      vmpi::mpi_complete_halo_swap();
      sim::calculate_spin_fields(post_comm_si, post_comm_ei);

      for (int atom = post_comm_si; atom < post_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k2_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k2_storage[atom].data(), half_dt);
      }

      for (int atom = pre_comm_si; atom < post_comm_ei; ++atom) writeback_spin(atom);

      //=====================================================================
      // K3 Stage (evaluate at t + dt/2)
      //=====================================================================
      vmpi::mpi_init_halo_swap();
      sim::calculate_spin_fields(pre_comm_si, pre_comm_ei);

      for (int atom = pre_comm_si; atom < pre_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k3_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k3_storage[atom].data(), dt);
      }

      vmpi::mpi_complete_halo_swap();
      sim::calculate_spin_fields(post_comm_si, post_comm_ei);

      for (int atom = post_comm_si; atom < post_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k3_storage[atom].data(), imaterial);
         predict_and_renorm(atom, k3_storage[atom].data(), dt);
      }

      for (int atom = pre_comm_si; atom < post_comm_ei; ++atom) writeback_spin(atom);

      //=====================================================================
      // K4 Stage (evaluate at t + dt) + final combine
      //
      // The final RK4 combination is folded into the per-atom K4 loop so
      // the comm-overlap pattern can do useful work on core atoms while
      // boundary halos are in flight.
      //=====================================================================
      vmpi::mpi_init_halo_swap();
      sim::calculate_spin_fields(pre_comm_si, pre_comm_ei);

      for (int atom = pre_comm_si; atom < pre_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k4_storage[atom].data(), imaterial);
         rk4_combine_and_renorm(atom, dt_over_6);
      }

      vmpi::mpi_complete_halo_swap();
      sim::calculate_spin_fields(post_comm_si, post_comm_ei);

      for (int atom = post_comm_si; atom < post_comm_ei; ++atom) {
         const int imaterial = atoms::type_array[atom];
         collect_H_FFT(atom, H);
         LL_FFT_method(y_pred_storage[atom].data(), H, k4_storage[atom].data(), imaterial);
         rk4_combine_and_renorm(atom, dt_over_6);
      }

      // Final write-back of all 9 components (core + boundary)
      for (int atom = pre_comm_si; atom < post_comm_ei; ++atom) writeback_full_state(atom);

      // Advance fine-grained noise time index
      quantum::increment_time();
   }

} // end of internal namespace
} // end of quantum namespace

#endif // MPICF
