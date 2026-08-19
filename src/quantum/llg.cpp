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
//   Top-level LLG dispatcher for the quantum thermostat module.
//
//   Selects the appropriate time stepper based on the chosen noise method
//   (FFT or HO) and the compilation mode (serial / MPI / CUDA).
//
//------------------------------------------------------------------------------

// Vampire headers
#include "quantum.hpp"

// Module headers
#include "internal.hpp"

#ifdef CUDA
#include "llg-HO-cuda.hpp"
#endif

namespace quantum{

   //------------------------------------------------------------------------
   // Public LLG dispatcher
   //------------------------------------------------------------------------
   void llg(){

      #ifdef MPICF
         // MPI parallel version
         if (internal::llg_method == internal::llg_fft) {
            internal::llg_FFT_mpi();
         }
         else if (internal::llg_method == internal::llg_ho) {
            internal::llg_HO_mpi();
         }
      #else
         #ifdef CUDA
            // CUDA accelerated version
            if (internal::llg_method == internal::llg_ho) {
               internal::cuda_ho::llg_HO_step();
            }
            else if (internal::llg_method == internal::llg_fft) {
               internal::llg_FFT();
            }
         #else
            // CPU serial version
            if (internal::llg_method == internal::llg_fft) {
               internal::llg_FFT();
            }
            else if (internal::llg_method == internal::llg_ho) {
               internal::llg_HO();
            }
         #endif
      #endif
   }

} // end of quantum namespace
