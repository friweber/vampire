//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2026. All rights reserved.
//
//------------------------------------------------------------------------------
//
//   FDT monitor (quantum:fdt-monitor, src/quantum/fdt_monitor.cpp) on the GPU:
//   gathers spin and external field of the probe atoms into one small device
//   buffer and copies it to the host. The external field of llg-heun holds the
//   thermal field (white or quantum bath) plus the applied fields; the applied
//   part is subtracted on the host (the monitor refuses dipole and FMR fields).
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"
#include "typedefs.hpp"

#ifdef CUDA

namespace cu = vcuda::internal;

namespace vcuda{

   namespace internal{
      namespace fdt{
         int* d_probe = 0;
         double* d_buf = 0;
         int n_probe = 0;

         __global__ void gather_kernel(const int* probe, const int np,
                                       const cu_real_t* sx, const cu_real_t* sy, const cu_real_t* sz,
                                       const cu_real_t* hx, const cu_real_t* hy, const cu_real_t* hz, double* buf){
            for(int p = blockIdx.x*blockDim.x + threadIdx.x; p < np; p += blockDim.x*gridDim.x){
               const int a = probe[p];
               buf[6*p]   = sx[a]; buf[6*p+1] = sy[a]; buf[6*p+2] = sz[a];
               buf[6*p+3] = hx[a]; buf[6*p+4] = hy[a]; buf[6*p+5] = hz[a];
            }
         }
      }
   }

   //---------------------------------------------------------------------------
   // spins s[3*np] and external fields h[3*np] of the probe atoms
   //---------------------------------------------------------------------------
   void fdt_gather(const std::vector<int>& probe, std::vector<double>& s, std::vector<double>& h){
      namespace f = cu::fdt;
      const int np = probe.size();
      if(f::d_probe == 0 || f::n_probe != np){
         if(f::d_probe){ cudaFree(f::d_probe); cudaFree(f::d_buf); }
         cudaMalloc((void**)&f::d_probe, np*sizeof(int));
         cudaMalloc((void**)&f::d_buf, 6*np*sizeof(double));
         cudaMemcpy(f::d_probe, probe.data(), np*sizeof(int), cudaMemcpyHostToDevice);
         f::n_probe = np;
      }
      f::gather_kernel <<< (np + 127)/128, 128 >>> (f::d_probe, np,
            cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
            cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field, f::d_buf);
      std::vector<double> buf(6*np);
      cudaMemcpy(buf.data(), f::d_buf, 6*np*sizeof(double), cudaMemcpyDeviceToHost);
      check_cuda_errors(__FILE__, __LINE__);
      for(int p = 0; p < np; p++)
         for(int c = 0; c < 3; c++){ s[3*p+c] = buf[6*p+c]; h[3*p+c] = buf[6*p+3+c]; }
   }

} // end of vcuda namespace

#endif
