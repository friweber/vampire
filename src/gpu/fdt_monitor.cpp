//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"

namespace gpu{

   //---------------------------------------------------------------------------
   // FDT monitor: spins and external fields of the probe atoms from the device
   //---------------------------------------------------------------------------
   void fdt_gather(const std::vector<int>& probe, std::vector<double>& s, std::vector<double>& h){
      #ifdef CUDA
         vcuda::fdt_gather(probe, s, h);
      #endif
      return;
   }

} // end of namespace gpu
