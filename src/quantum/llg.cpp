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

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"
#include "quantum.hpp"

// quantum module headers
#include "internal.hpp"
#include "llg_atom.hpp"

namespace quantum{

   //---------------------------------------------------------------------------
   // Function to integrate one open-system LLG step. Selects the stepper for
   // the build (serial or MPI) and the noise generator. A GPU stepper would
   // be dispatched here, ahead of the CPU choice, following gpu::llg_heun().
   //---------------------------------------------------------------------------
   void llg(){

      using namespace internal;

      #ifdef MPICF
         if(noise_generator == on_the_fly) llg_ho_mpi();
         else llg_fft_mpi();
      #else
         if(noise_generator == on_the_fly) llg_ho_serial();
         else llg_fft_serial();
      #endif

      return;

   }

   namespace internal{

      //---------------------------------------------------------------------------
      // Function to allocate the auxiliary oscillator and the RK4 storage
      //---------------------------------------------------------------------------
      void allocate_thermostat(const int n_atoms){

         q_x.assign(n_atoms, 0.0);
         q_y.assign(n_atoms, 0.0);
         q_z.assign(n_atoms, 0.0);
         p_x.assign(n_atoms, 0.0);
         p_y.assign(n_atoms, 0.0);
         p_z.assign(n_atoms, 0.0);

         const size_t n = 9*static_cast<size_t>(n_atoms);
         k1.assign(n, 0.0);
         k2.assign(n, 0.0);
         k3.assign(n, 0.0);
         k4.assign(n, 0.0);
         y_pred.assign(n, 0.0);
         y_in.assign(n, 0.0);

         return;

      }

      //---------------------------------------------------------------------------
      // Function to copy the state (S, q, p) at the start of the step
      //---------------------------------------------------------------------------
      void rk4_save(const int start_index, const int end_index){

         for(int atom = start_index; atom < end_index; atom++){
            double* y = &y_in[9*atom];
            y[0] = atoms::x_spin_array[atom];
            y[1] = atoms::y_spin_array[atom];
            y[2] = atoms::z_spin_array[atom];
            y[3] = q_x[atom];
            y[4] = q_y[atom];
            y[5] = q_z[atom];
            y[6] = p_x[atom];
            y[7] = p_y[atom];
            y[8] = p_z[atom];
         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to evaluate one RK4 stage for a range of atoms. Stage 1 uses
      // the saved state, stages 2 to 4 the predictor; stages 1 to 3 also form
      // the next predictor. ho selects where this step's bath sample enters:
      // the oscillator (on-the-fly noise) or the field (pre-generated noise).
      //---------------------------------------------------------------------------
      void rk4_stage(const int start_index, const int end_index, const int stage, const bool ho){

         const double dt = mp::dt;
         const double coeff = (stage == 3) ? dt : 0.5*dt;

         std::vector<double>& k = (stage == 1) ? k1 : (stage == 2) ? k2 : (stage == 3) ? k3 : k4;
         const std::vector<double>& y_src = (stage == 1) ? y_in : y_pred;
         const bath_t& b = spin_bath;

         double H[3];
         double noise[3];

         for(int atom = start_index; atom < end_index; atom++){

            const int m = atoms::type_array[atom];
            const double* y = &y_src[9*atom];
            double* dydt = &k[9*atom];

            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom];
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom];
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom];

            if(ho){
               noise[0] = b.x[atom];
               noise[1] = b.y[atom];
               noise[2] = b.z[atom];
               ho_equation_of_motion(y, H, noise, material_A[m], material_gamma[m], material_omega0[m], dydt);
            }
            else{
               H[0] += b.x[atom];
               H[1] += b.y[atom];
               H[2] += b.z[atom];
               fft_equation_of_motion(y, H, material_A[m], material_gamma[m], material_omega0[m], dydt);
            }

            if(stage < 4) rk4_predict(&y_in[9*atom], dydt, coeff, &y_pred[9*atom]);

         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to publish the predicted spins so the fields can be recomputed
      //---------------------------------------------------------------------------
      void rk4_writeback_spin(const int start_index, const int end_index){

         for(int atom = start_index; atom < end_index; atom++){
            atoms::x_spin_array[atom] = y_pred[9*atom + 0];
            atoms::y_spin_array[atom] = y_pred[9*atom + 1];
            atoms::z_spin_array[atom] = y_pred[9*atom + 2];
         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to combine the four stages and commit the new state
      //---------------------------------------------------------------------------
      void rk4_finish(const int start_index, const int end_index){

         const double dt_over_6 = mp::dt/6.0;

         for(int atom = start_index; atom < end_index; atom++){
            double* y = &y_pred[9*atom];
            rk4_combine(&y_in[9*atom], &k1[9*atom], &k2[9*atom], &k3[9*atom], &k4[9*atom], dt_over_6, y);
            atoms::x_spin_array[atom] = y[0];
            atoms::y_spin_array[atom] = y[1];
            atoms::z_spin_array[atom] = y[2];
            q_x[atom] = y[3];
            q_y[atom] = y[4];
            q_z[atom] = y[5];
            p_x[atom] = y[6];
            p_y[atom] = y[7];
            p_z[atom] = y[8];
         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to export the auxiliary oscillator of one atom, the field the
      // spin sees from the bath
      //---------------------------------------------------------------------------
      void export_thermostat_sample(){

         if(export_noise) export_sample(q_x[export_atom], q_y[export_atom], q_z[export_atom]);

         return;

      }

   } // end of internal namespace

} // end of quantum namespace
