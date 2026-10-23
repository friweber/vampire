//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//------------------------------------------------------------------------------

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"

namespace gpu{

   //---------------------------------------------------------------------------
   // Open-system LLG (quantum module, llg-quantum) on the GPU
   //---------------------------------------------------------------------------
   void quantum_llg(){
      #ifdef CUDA
         vcuda::quantum_llg();
      #endif
      return;
   }

} // end of namespace gpu
