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

#ifdef MPICF

// C++ standard library headers

// Vampire headers
#include "quantum.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      //---------------------------------------------------------------------------
      // One RK4 stage across core and boundary atoms. The halo swap for the
      // boundary atoms overlaps with the core evaluation; the external field,
      // which does not change within the step, is only computed at stage 1.
      //---------------------------------------------------------------------------
      static void mpi_stage(const int stage, const bool ho){

         const int core = vmpi::num_core_atoms;
         const int n = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         vmpi::mpi_init_halo_swap();
         sim::calculate_spin_fields(0, core);
         if(stage == 1) sim::calculate_external_fields(0, n);
         rk4_stage(0, core, stage, ho);

         vmpi::mpi_complete_halo_swap();
         sim::calculate_spin_fields(core, n);
         if(stage == 1) sim::calculate_external_fields(core, n);
         rk4_stage(core, n, stage, ho);

         if(stage < 4) rk4_writeback_spin(0, n);

         return;

      }

      //---------------------------------------------------------------------------
      // MPI RK4 step with on-the-fly noise. The bath sample is drawn once for
      // all local atoms and drives the auxiliary oscillator in every stage.
      //---------------------------------------------------------------------------
      void llg_ho_mpi(){

         const int n = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         vmpi::mpi_init_halo_swap();
         vmpi::mpi_complete_halo_swap();

         rk4_save(0, n);
         draw(spin_bath);

         mpi_stage(1, true);
         mpi_stage(2, true);
         mpi_stage(3, true);
         mpi_stage(4, true);

         rk4_finish(0, n);

         export_thermostat_sample();

         vmpi::barrier();

         return;

      }

      //---------------------------------------------------------------------------
      // MPI RK4 step with pre-generated noise. The sample is already shaped
      // by the Lorentzian and is added to the field in every stage.
      //---------------------------------------------------------------------------
      void llg_fft_mpi(){

         const int n = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         vmpi::mpi_init_halo_swap();
         vmpi::mpi_complete_halo_swap();

         rk4_save(0, n);
         draw(spin_bath);

         mpi_stage(1, false);
         mpi_stage(2, false);
         mpi_stage(3, false);
         mpi_stage(4, false);

         rk4_finish(0, n);

         export_thermostat_sample();

         vmpi::barrier();

         return;

      }

   } // end of internal namespace

} // end of quantum namespace

#endif // MPICF
