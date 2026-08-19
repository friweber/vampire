//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
// MPI-parallel Heun LLG integrator with quantum colored noise.
//
// Mirrors LLGHeun-mpi.cpp structure with core/boundary split and halo swaps.
// Quantum spin noise is injected after each calculate_external_fields() call
// (once for core Euler, once for boundary Euler, once for core Heun, once for
// boundary Heun), using the same pre-drawn values per step — correct treatment.
//
//-----------------------------------------------------------------------------
#ifdef MPICF
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "quantum.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "vmpi.hpp"

// Module internal headers
#include "../quantum/internal.hpp"        // quantum::internal::heun_noise_kind
#include "../spinlattice/internal.hpp"    // sld::internal::export_noise

#include <cmath>

int calculate_spin_fields(const int, const int);
int calculate_external_fields(const int, const int);

namespace sim{

int LLG_Heun_quantum_mpi(){

   using namespace LLG_arrays;

   // Check for initialisation of LLG integration arrays
   if(LLG_set==false) sim::LLGinit();

   // -------------------------------------------------------------------------
   // First-call initialisation of quantum noise bath
   // -------------------------------------------------------------------------
   static bool q_heun_mpi_initialized = false;
   if (!q_heun_mpi_initialized) {
      q_heun_mpi_initialized = true;

      sim::hamiltonian_simulation_flags[3] = 0;

      quantum::sld_noise::kind_t kind = quantum::internal::heun_noise_kind;
      if (sld::enabled && sld::internal::quantum_noise_type != sld::internal::sld_classical) {
         switch (sld::internal::quantum_noise_type) {
            case sld::internal::sld_quantum:
               kind = quantum::sld_noise::kind_t::quantum;               break;
            case sld::internal::sld_quantum_no_zero:
               kind = quantum::sld_noise::kind_t::quantum_no_zero;       break;
            case sld::internal::sld_quantum_fft:
               kind = quantum::sld_noise::kind_t::quantum_fft;           break;
            case sld::internal::sld_quantum_no_zero_fft:
               kind = quantum::sld_noise::kind_t::quantum_no_zero_fft;   break;
            default:
               kind = quantum::sld_noise::kind_t::quantum;               break;
         }
      }

      const bool use_fft = (kind == quantum::sld_noise::kind_t::quantum_fft ||
                            kind == quantum::sld_noise::kind_t::quantum_no_zero_fft);
      const uint64_t n_fine = use_fft
         ? (sim::equilibration_time + sim::total_time)
         : 0;

      const int num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      std::vector<quantum::sld_noise::material_params> mats(mp::num_materials);
      for (int m = 0; m < mp::num_materials; ++m) {
         mats[m].mass       = 0.0;
         mats[m].damp_lat   = 0.0;
         mats[m].V0         = 0.0;
         mats[m].H_th_sigma = mp::material[m].H_th_sigma;
      }

      quantum::sld_noise::initialize(kind, num_atoms, mats, n_fine);

      if (sld::internal::export_noise && vmpi::my_rank == 0) {
         quantum::sld_noise::enable_spin_noise_export(
            sld::internal::export_noise_filename,
            sld::internal::export_noise_atom);
      }
   }

   // -------------------------------------------------------------------------
   // Per-step integration
   // -------------------------------------------------------------------------
   const int pre_comm_si  = 0;
   const int pre_comm_ei  = vmpi::num_core_atoms;
   const int post_comm_si = vmpi::num_core_atoms;
   const int post_comm_ei = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
   const int num_atoms    = post_comm_ei;

   double xyz[3];
   double S_new[3];
   double mod_S;

   // Draw new quantum noise values for this step
   quantum::sld_noise::generate(num_atoms);

   // Initiate halo swap
   vmpi::mpi_init_halo_swap();

   // Store initial spin positions (all)
   for(int atom = pre_comm_si; atom < post_comm_ei; atom++){
      x_initial_spin_array[atom] = atoms::x_spin_array[atom];
      y_initial_spin_array[atom] = atoms::y_spin_array[atom];
      z_initial_spin_array[atom] = atoms::z_spin_array[atom];
   }

   // Calculate fields (core)
   calculate_spin_fields(pre_comm_si, pre_comm_ei);
   calculate_external_fields(pre_comm_si, pre_comm_ei);

   // Inject quantum spin noise (core)
   for(int atom = pre_comm_si; atom < pre_comm_ei; atom++){
      atoms::x_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 0);
      atoms::y_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 1);
      atoms::z_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 2);
   }

   // Calculate Euler Step (core)
   for(int atom = pre_comm_si; atom < pre_comm_ei; atom++){

      const int imaterial = atoms::type_array[atom];
      const double one_oneplusalpha_sq   = material_parameters::material[imaterial].one_oneplusalpha_sq;
      const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

      const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
      const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

      xyz[0] = (one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
      xyz[1] = (one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
      xyz[2] = (one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

      x_euler_array[atom] = xyz[0];
      y_euler_array[atom] = xyz[1];
      z_euler_array[atom] = xyz[2];

      S_new[0] = S[0] + xyz[0]*material_parameters::dt;
      S_new[1] = S[1] + xyz[1]*material_parameters::dt;
      S_new[2] = S[2] + xyz[2]*material_parameters::dt;

      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

      x_spin_storage_array[atom] = S_new[0]*mod_S;
      y_spin_storage_array[atom] = S_new[1]*mod_S;
      z_spin_storage_array[atom] = S_new[2]*mod_S;
   }

   // Complete halo swap
   vmpi::mpi_complete_halo_swap();

   // Calculate fields (boundary)
   calculate_spin_fields(post_comm_si, post_comm_ei);
   calculate_external_fields(post_comm_si, post_comm_ei);

   // Inject quantum spin noise (boundary)
   for(int atom = post_comm_si; atom < post_comm_ei; atom++){
      atoms::x_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 0);
      atoms::y_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 1);
      atoms::z_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 2);
   }

   // Calculate Euler Step (boundary)
   for(int atom = post_comm_si; atom < post_comm_ei; atom++){

      const int imaterial = atoms::type_array[atom];
      const double one_oneplusalpha_sq   = material_parameters::material[imaterial].one_oneplusalpha_sq;
      const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

      const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
      const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

      xyz[0] = (one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
      xyz[1] = (one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
      xyz[2] = (one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

      x_euler_array[atom] = xyz[0];
      y_euler_array[atom] = xyz[1];
      z_euler_array[atom] = xyz[2];

      S_new[0] = S[0] + xyz[0]*material_parameters::dt;
      S_new[1] = S[1] + xyz[1]*material_parameters::dt;
      S_new[2] = S[2] + xyz[2]*material_parameters::dt;

      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

      x_spin_storage_array[atom] = S_new[0]*mod_S;
      y_spin_storage_array[atom] = S_new[1]*mod_S;
      z_spin_storage_array[atom] = S_new[2]*mod_S;
   }

   // Copy new spins to spin array (all)
   for(int atom = pre_comm_si; atom < post_comm_ei; atom++){
      atoms::x_spin_array[atom] = x_spin_storage_array[atom];
      atoms::y_spin_array[atom] = y_spin_storage_array[atom];
      atoms::z_spin_array[atom] = z_spin_storage_array[atom];
   }

   // Initiate second halo swap
   vmpi::mpi_init_halo_swap();

   // Recalculate spin dependent fields (core)
   calculate_spin_fields(pre_comm_si, pre_comm_ei);
   calculate_external_fields(pre_comm_si, pre_comm_ei);

   // Inject quantum spin noise for Heun correction (core)
   for(int atom = pre_comm_si; atom < pre_comm_ei; atom++){
      atoms::x_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 0);
      atoms::y_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 1);
      atoms::z_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 2);
   }

   // Calculate Heun Gradients (core)
   for(int atom = pre_comm_si; atom < pre_comm_ei; atom++){

      const int imaterial = atoms::type_array[atom];
      const double one_oneplusalpha_sq   = material_parameters::material[imaterial].one_oneplusalpha_sq;
      const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

      const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
      const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

      xyz[0] = (one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
      xyz[1] = (one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
      xyz[2] = (one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

      x_heun_array[atom] = xyz[0];
      y_heun_array[atom] = xyz[1];
      z_heun_array[atom] = xyz[2];
   }

   // Complete second halo swap
   vmpi::mpi_complete_halo_swap();

   // Recalculate spin dependent fields (boundary)
   calculate_spin_fields(post_comm_si, post_comm_ei);
   calculate_external_fields(post_comm_si, post_comm_ei);

   // Inject quantum spin noise for Heun correction (boundary)
   for(int atom = post_comm_si; atom < post_comm_ei; atom++){
      atoms::x_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 0);
      atoms::y_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 1);
      atoms::z_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 2);
   }

   // Calculate Heun Gradients (boundary)
   for(int atom = post_comm_si; atom < post_comm_ei; atom++){

      const int imaterial = atoms::type_array[atom];
      const double one_oneplusalpha_sq   = material_parameters::material[imaterial].one_oneplusalpha_sq;
      const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

      const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
      const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

      xyz[0] = (one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
      xyz[1] = (one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
      xyz[2] = (one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

      x_heun_array[atom] = xyz[0];
      y_heun_array[atom] = xyz[1];
      z_heun_array[atom] = xyz[2];
   }

   // Calculate Heun Step (all)
   for(int atom = pre_comm_si; atom < post_comm_ei; atom++){
      S_new[0] = x_initial_spin_array[atom] + material_parameters::half_dt*(x_euler_array[atom] + x_heun_array[atom]);
      S_new[1] = y_initial_spin_array[atom] + material_parameters::half_dt*(y_euler_array[atom] + y_heun_array[atom]);
      S_new[2] = z_initial_spin_array[atom] + material_parameters::half_dt*(z_euler_array[atom] + z_heun_array[atom]);

      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

      atoms::x_spin_array[atom] = S_new[0]*mod_S;
      atoms::y_spin_array[atom] = S_new[1]*mod_S;
      atoms::z_spin_array[atom] = S_new[2]*mod_S;
   }

   // Optional noise export (rank 0 only)
   if (sld::internal::export_noise && vmpi::my_rank == 0) {
      quantum::sld_noise::export_spin_noise_step(sim::time * mp::dt_SI);
   }

   // Swap timers compute -> wait
   vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

   vmpi::barrier();

   // Swap timers wait -> compute
   vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

   return EXIT_SUCCESS;
}

} // end of namespace sim
#endif
