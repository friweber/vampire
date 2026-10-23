//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2026. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//
//   Quantum module on the GPU
//
//   Coloured spin bath for llg-heun: the on-the-fly (HO) generators of the
//   quantum module -- Matsubara Ornstein-Uhlenbeck modes for "quantum",
//   log-spaced OU modes for "quantum-no-zero" -- advanced once per step on
//   the device, one thread per atom. The sample replaces the white thermal
//   field and is seen by predictor and corrector alike.
//
//   The host keeps the bath coefficients (rebuilt when sim::temperature moves,
//   lambda_k rescaled with T for the log bath) and uploads them to constant
//   memory only when they change. Random numbers: one Philox stream per thread.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"
#include "typedefs.hpp"

#ifdef CUDA

#include <curand_kernel.h>

namespace cu = vcuda::internal;

namespace vcuda{

   namespace internal{

      namespace quantum_gpu{

         constexpr int max_modes = 256;
         constexpr int max_materials = 64;

         // bath coefficients (see gpu_view_t)
         __constant__ double c_c0[max_modes];
         __constant__ double c_c1[max_modes];
         __constant__ double c_c2[max_modes];
         __constant__ double c_c3[max_modes];
         __constant__ double c_amplitude[max_materials];

         bool initialised = false;
         int n_atoms = 0;
         int n_modes = 0;
         int spectrum = 0;
         double uploaded_T = -1.0;
         const double* uploaded_amplitude = 0;

         curandStatePhilox4_32_10_t* d_rng = 0;
         double* d_s = 0;      // OU mode states, [component][mode][atom]
         double* d_noise = 0;  // this step's sample, [component][atom]

         //---------------------------------------------------------------------
         // Kernels
         //---------------------------------------------------------------------
         __global__ void init_rng_kernel(curandStatePhilox4_32_10_t* states, unsigned long long seed, int n){
            const int tid = blockIdx.x*blockDim.x + threadIdx.x;
            if(tid < n) curand_init(seed, tid, 0, &states[tid]);
         }

         // Gaussian draws: single precision, four per Philox call. FP64 is slow on consumer
         // GPUs (1/64 rate on GeForce) and the draws need no double precision; the OU
         // recursion and the sums stay in double.
         struct normal_source_t{
            curandStatePhilox4_32_10_t* rng;
            float4 buf;
            int left;
            __device__ normal_source_t(curandStatePhilox4_32_10_t* r) : rng(r), left(0) {}
            __device__ double next(){
               if(left == 0){ buf = curand_normal4(rng); left = 4; }
               left--;
               return static_cast<double>(left == 3 ? buf.x : left == 2 ? buf.y : left == 1 ? buf.z : buf.w);
            }
         };

         // One bath step for every atom. spectrum 1: Matsubara OU (+ white zero mode),
         // spectrum 2: log-spaced OU. Same update as ou_update_atom / log_bath_update_atom.
         __global__ void bath_kernel(double* s, double* out, const int* material,
                                     curandStatePhilox4_32_10_t* states,
                                     const int n_atoms, const int n_modes, const int spectrum,
                                     const double white_amp){
            const int tid = blockIdx.x*blockDim.x + threadIdx.x;
            curandStatePhilox4_32_10_t rng = states[tid];
            normal_source_t gauss(&rng);
            const size_t N = n_atoms;
            for(int atom = tid; atom < n_atoms; atom += blockDim.x*gridDim.x){
               const double amplitude = c_amplitude[material[atom]];
               for(int c = 0; c < 3; c++){
                  double* sc = s + static_cast<size_t>(c)*n_modes*N;
                  double sum = 0.0;
                  if(spectrum == 1){
                     for(int k = 0; k < n_modes; k++){
                        const double xi = gauss.next();
                        double sk = sc[k*N + atom];
                        sum += c_c2[k]*sk + c_c3[k]*xi;        // drift*s(old) + noise_amp*xi
                        sk = c_c0[k]*sk + c_c1[k]*xi;          // exact OU update
                        sc[k*N + atom] = sk;
                     }
                     sum += white_amp*gauss.next();
                  }
                  else{
                     for(int k = 0; k < n_modes; k++){
                        const double xi = gauss.next();
                        const double sk = c_c0[k]*sc[k*N + atom] + c_c1[k]*xi;
                        sc[k*N + atom] = sk;
                        sum += c_c2[k]*sk;                      // amp * s(new)
                     }
                  }
                  out[c*N + atom] = amplitude*sum;
               }
            }
            states[tid] = rng;
         }

         // Copy the bath sample into the thermal (external) field arrays
         __global__ void bath_to_field_kernel(const double* noise, cu_real_t* hx, cu_real_t* hy, cu_real_t* hz, const int n_atoms){
            const size_t N = n_atoms;
            for(int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < n_atoms; atom += blockDim.x*gridDim.x){
               hx[atom] = noise[atom];
               hy[atom] = noise[N + atom];
               hz[atom] = noise[2*N + atom];
            }
         }

         //---------------------------------------------------------------------
         // Host functions
         //---------------------------------------------------------------------
         void initialise(){

            quantum::gpu_view_t v;
            quantum::gpu_view(v);

            n_atoms = ::atoms::num_atoms;
            n_modes = v.n_modes;
            spectrum = v.spectrum;

            if(n_modes > max_modes) err::fatal("Error - quantum:bath-modes exceeds the GPU limit of 256.");
            if(v.n_materials > max_materials) err::fatal("Error - the quantum module on the GPU supports at most 64 materials.");

            const int n_threads = cu::grid_size*cu::block_size;
            cudaMalloc((void**)&d_rng, n_threads*sizeof(curandStatePhilox4_32_10_t));
            init_rng_kernel <<< cu::grid_size, cu::block_size >>> (d_rng, 0x9e3779b97f4a7c15ULL ^ static_cast<unsigned long long>(::mtrandom::integration_seed), n_threads);

            const size_t N = n_atoms;
            cudaMalloc((void**)&d_noise, 3*N*sizeof(double));
            cudaMemset(d_noise, 0, 3*N*sizeof(double));
            if(spectrum != 0){
               cudaMalloc((void**)&d_s, 3*N*n_modes*sizeof(double));
               cudaMemset(d_s, 0, 3*N*n_modes*sizeof(double));   // as the CPU bath: modes start from zero
            }
            check_cuda_errors(__FILE__, __LINE__);

            zlog << zTs() << "Quantum module on the GPU: " << "coloured bath for llg-heun"
                 << ", spectrum " << (spectrum == 0 ? "classical" : spectrum == 1 ? "quantum (Matsubara OU)" : "quantum-no-zero (log OU)")
                 << ", " << n_modes << " modes, " << n_atoms << " atoms" << std::endl;

            initialised = true;
            return;
         }

         // Upload coefficients when the temperature or the phase changed; return the view
         void prepare(quantum::gpu_view_t& v){
            if(!initialised) initialise();
            quantum::gpu_view(v);
            if(v.spectrum != 0 && v.T_scaled != uploaded_T){
               cudaMemcpyToSymbol(c_c0, v.c0, n_modes*sizeof(double));
               cudaMemcpyToSymbol(c_c1, v.c1, n_modes*sizeof(double));
               cudaMemcpyToSymbol(c_c2, v.c2, n_modes*sizeof(double));
               if(v.c3) cudaMemcpyToSymbol(c_c3, v.c3, n_modes*sizeof(double));
               uploaded_T = v.T_scaled;
            }
            if(v.amplitude != uploaded_amplitude){
               cudaMemcpyToSymbol(c_amplitude, v.amplitude, v.n_materials*sizeof(double));
               uploaded_amplitude = v.amplitude;
            }
            return;
         }

         void draw(const quantum::gpu_view_t& v){
            bath_kernel <<< cu::grid_size, cu::block_size >>> (d_s, d_noise, cu::atoms::d_materials, d_rng,
                                                                n_atoms, n_modes, v.spectrum, v.white_amp);
            check_cuda_errors(__FILE__, __LINE__);
         }

      } // end of quantum_gpu namespace

      //------------------------------------------------------------------------
      // Coloured thermal field for llg-heun (called from update_external_fields)
      //------------------------------------------------------------------------
      void quantum_bath_fields(){
         namespace qg = quantum_gpu;
         quantum::gpu_view_t v;
         qg::prepare(v);
         qg::draw(v);
         qg::bath_to_field_kernel <<< cu::grid_size, cu::block_size >>> (qg::d_noise,
               cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field, qg::n_atoms);
         check_cuda_errors(__FILE__, __LINE__);
      }

      void quantum_finalize(){
         namespace qg = quantum_gpu;
         if(!qg::initialised) return;
         cudaFree(qg::d_rng); cudaFree(qg::d_noise);
         if(qg::d_s) cudaFree(qg::d_s);
         qg::initialised = false;
      }

   } // end of internal namespace

} // end of vcuda namespace

#endif
