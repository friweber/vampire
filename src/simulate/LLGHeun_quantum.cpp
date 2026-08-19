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
// Heun LLG integrator with quantum colored noise.
//
// Noise amplitude is derived from H_th_sigma only (no quantum-lorentzian-*
// parameters needed). Noise kind is selected via spin-lattice:noise-type.
//
// The standard thermal field is disabled and replaced by quantum spin noise
// injected into the external field array after each calculate_external_fields()
// call (once for the Euler step, once for the Heun correction, using the same
// pre-drawn values — correct Heun treatment of multiplicative noise).
//
//-----------------------------------------------------------------------------

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "sim.hpp"
#include "sld.hpp"

// Module internal headers
#include "../quantum/internal.hpp"        // quantum::internal::heun_noise_kind
#include "../spinlattice/internal.hpp"    // sld::internal::export_noise

namespace sim{

int LLG_Heun_quantum(){

   using namespace LLG_arrays;

   // Check for initialisation of LLG integration arrays
   if(LLG_set==false) sim::LLGinit();

   // -------------------------------------------------------------------------
   // First-call initialisation of quantum noise bath
   // -------------------------------------------------------------------------
   static bool q_heun_initialized = false;
   if (!q_heun_initialized) {
      q_heun_initialized = true;

      // Disable classical thermal field — quantum noise provides it instead
      sim::hamiltonian_simulation_flags[3] = 0;

      // Noise kind: quantum:heun-noise-type (ASD path, no SLD overhead).
      // Falls back to spin-lattice:noise-type for FFT variants when SLD is active.
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

      // Build per-material params. mass=0 disables phonon bath (ASD, no lattice).
      std::vector<quantum::sld_noise::material_params> mats(mp::num_materials);
      for (int m = 0; m < mp::num_materials; ++m) {
         mats[m].mass       = 0.0;
         mats[m].damp_lat   = 0.0;
         mats[m].V0         = 0.0;
         mats[m].H_th_sigma = mp::material[m].H_th_sigma;
      }

      quantum::sld_noise::initialize(kind, atoms::num_atoms, mats, n_fine);

      if (sld::internal::export_noise) {
         quantum::sld_noise::enable_spin_noise_export(
            sld::internal::export_noise_filename,
            sld::internal::export_noise_atom);
      }
   }

   // -------------------------------------------------------------------------
   // Per-step integration
   // -------------------------------------------------------------------------
   const int num_atoms = atoms::num_atoms;
   double xyz[3];
   double S_new[3];
   double mod_S;

   // Draw new quantum noise values for this step
   quantum::sld_noise::generate(num_atoms);

   // Store initial spin positions
   for(int atom = 0; atom < num_atoms; atom++){
      x_initial_spin_array[atom] = atoms::x_spin_array[atom];
      y_initial_spin_array[atom] = atoms::y_spin_array[atom];
      z_initial_spin_array[atom] = atoms::z_spin_array[atom];
   }

   // Calculate fields (zeros external field array first internally)
   calculate_spin_fields(0, num_atoms);
   calculate_external_fields(0, num_atoms);

   // Inject quantum spin noise into external field array
   for(int atom = 0; atom < num_atoms; atom++){
      atoms::x_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 0);
      atoms::y_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 1);
      atoms::z_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 2);
   }

   // Calculate Euler Step
   for(int atom = 0; atom < num_atoms; atom++){

      const int imaterial = atoms::type_array[atom];
      const double one_oneplusalpha_sq   = mp::material[imaterial].one_oneplusalpha_sq;
      const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;

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

      S_new[0] = S[0] + xyz[0]*mp::dt;
      S_new[1] = S[1] + xyz[1]*mp::dt;
      S_new[2] = S[2] + xyz[2]*mp::dt;

      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

      x_spin_storage_array[atom] = S_new[0]*mod_S;
      y_spin_storage_array[atom] = S_new[1]*mod_S;
      z_spin_storage_array[atom] = S_new[2]*mod_S;
   }

   // Copy Euler spins
   for(int atom = 0; atom < num_atoms; atom++){
      atoms::x_spin_array[atom] = x_spin_storage_array[atom];
      atoms::y_spin_array[atom] = y_spin_storage_array[atom];
      atoms::z_spin_array[atom] = z_spin_storage_array[atom];
   }

   // Recalculate spin dependent fields; re-inject same quantum noise (Heun)
   calculate_spin_fields(0, num_atoms);
   calculate_external_fields(0, num_atoms);

   for(int atom = 0; atom < num_atoms; atom++){
      atoms::x_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 0);
      atoms::y_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 1);
      atoms::z_total_external_field_array[atom] += quantum::sld_noise::spin(atom, 2);
   }

   // Calculate Heun Gradients
   for(int atom = 0; atom < num_atoms; atom++){

      const int imaterial = atoms::type_array[atom];
      const double one_oneplusalpha_sq   = mp::material[imaterial].one_oneplusalpha_sq;
      const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;

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

   // Calculate Heun Step
   for(int atom = 0; atom < num_atoms; atom++){
      S_new[0] = x_initial_spin_array[atom] + mp::half_dt*(x_euler_array[atom] + x_heun_array[atom]);
      S_new[1] = y_initial_spin_array[atom] + mp::half_dt*(y_euler_array[atom] + y_heun_array[atom]);
      S_new[2] = z_initial_spin_array[atom] + mp::half_dt*(z_euler_array[atom] + z_heun_array[atom]);

      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

      atoms::x_spin_array[atom] = S_new[0]*mod_S;
      atoms::y_spin_array[atom] = S_new[1]*mod_S;
      atoms::z_spin_array[atom] = S_new[2]*mod_S;
   }

   // Optional noise export for this step
   if (sld::internal::export_noise) {
      quantum::sld_noise::export_spin_noise_step(sim::time * mp::dt_SI);
   }

   return EXIT_SUCCESS;
}

} // end of namespace sim
