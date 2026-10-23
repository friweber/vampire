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
#include "quantum.hpp"
#include "sim.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      //---------------------------------------------------------------------------
      // Serial RK4 step with on-the-fly noise. The bath sample is drawn once
      // and drives the auxiliary oscillator in every stage.
      //---------------------------------------------------------------------------
      void llg_ho_serial(){

         const int n = spin_bath.n_sites;

         rk4_save(0, n);
         sim::calculate_spin_fields(0, n);
         sim::calculate_external_fields(0, n);
         draw(spin_bath);

         rk4_stage(0, n, 1, true);
         rk4_writeback_spin(0, n);
         sim::calculate_spin_fields(0, n);

         rk4_stage(0, n, 2, true);
         rk4_writeback_spin(0, n);
         sim::calculate_spin_fields(0, n);

         rk4_stage(0, n, 3, true);
         rk4_writeback_spin(0, n);
         sim::calculate_spin_fields(0, n);

         rk4_stage(0, n, 4, true);
         rk4_finish(0, n);

         export_thermostat_sample();

         return;

      }

      //---------------------------------------------------------------------------
      // Serial RK4 step with pre-generated noise. The sample is already
      // shaped by the Lorentzian and is added to the field in every stage.
      //---------------------------------------------------------------------------
      void llg_fft_serial(){

         const int n = spin_bath.n_sites;

         rk4_save(0, n);
         sim::calculate_spin_fields(0, n);
         sim::calculate_external_fields(0, n);
         draw(spin_bath);

         rk4_stage(0, n, 1, false);
         rk4_writeback_spin(0, n);
         sim::calculate_spin_fields(0, n);

         rk4_stage(0, n, 2, false);
         rk4_writeback_spin(0, n);
         sim::calculate_spin_fields(0, n);

         rk4_stage(0, n, 3, false);
         rk4_writeback_spin(0, n);
         sim::calculate_spin_fields(0, n);

         rk4_stage(0, n, 4, false);
         rk4_finish(0, n);

         export_thermostat_sample();

         return;

      }

   } // end of internal namespace

} // end of quantum namespace
