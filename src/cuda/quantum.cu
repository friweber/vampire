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
//   1. Coloured spin bath for llg-heun: the on-the-fly (HO) generators of the
//      quantum module -- Matsubara Ornstein-Uhlenbeck modes for "quantum",
//      log-spaced OU modes for "quantum-no-zero" -- advanced once per step on
//      the device, one thread per atom. The sample replaces the white thermal
//      field and is seen by predictor and corrector alike.
//   2. Open-system LLG (llg-quantum): spin + Lorentzian auxiliary oscillator,
//      classical RK4, with the bath noise driving the oscillator. The per-atom
//      arithmetic is the module's own (src/quantum/llg_atom.hpp), so CPU and
//      GPU integrate the same equations. Combined thermostat (some alpha_e >
//      0): a second bath of the same spectrum, the Ohmic (electron) bath,
//      adds its field and a Gilbert damping alpha_e to the spin equation.
//
//   3. Pre-generated (FFT) noise, for both: llg-heun uses the sample as the
//      thermal field, the open-system LLG adds it to the spin field (the
//      Lorentzian is in its spectrum) and integrates the oscillator without
//      noise (fft_equation_of_motion). Two generators:
//      - GPU (default): the overlap-save windows of src/quantum/bath_fft.cpp
//        are generated on the device with cuFFT (white noise from the Philox
//        streams, batched D2Z, filter, Z2D) and interpolated per step by a
//        kernel; the window lives in device memory, the host does nothing.
//      - host (quantum:pre-generated-device = host): the CPU generator draws
//        the sample and it is uploaded once per step; the host RNG then draws
//        exactly what a CPU run with the same seed draws.
//
//   The host keeps the bath coefficients (rebuilt when sim::temperature moves,
//   lambda_k rescaled with T for the log bath) and uploads them to constant
//   memory only when they change. Random numbers: one Philox stream per thread.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <sstream>
#include <string>
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

// quantum module per-atom arithmetic (host and device)
#include "../quantum/llg_atom.hpp"

#ifdef CUDA

#include <curand_kernel.h>
#include <cufft.h>

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
         __constant__ double c_white[max_materials];
         __constant__ double c_A[max_materials];
         __constant__ double c_gamma[max_materials];
         __constant__ double c_omega0[max_materials];
         __constant__ double c_amplitude_e[max_materials]; // combined thermostat: Ohmic bath
         __constant__ double c_white_e[max_materials];
         __constant__ double c_alpha_e[max_materials];

         bool initialised = false;
         int n_atoms = 0;
         int n_modes = 0;
         int spectrum = 0;
         int mode = 0;
         int ohmic = 0;
         int pregen = 0;       // pre-generated (FFT) noise from the host
         double uploaded_T = -1.0;
         const double* uploaded_amplitude = 0;

         curandStatePhilox4_32_10_t* d_rng = 0;
         double* d_s = 0;      // OU mode states, [component][mode][atom]
         double* d_noise = 0;  // this step's sample, [component][atom]
         double* d_q = 0;      // auxiliary oscillator position, [component][atom]
         double* d_p = 0;      // auxiliary oscillator momentum
         double* d_yin = 0;    // RK4 state at the start of the step, [9][atom]
         double* d_ypred = 0;  // RK4 predictor
         double* d_k = 0;      // RK4 slopes, [4][9][atom]
         double* d_s_e = 0;    // combined thermostat: Ohmic bath mode states, as d_s
         double* d_noise_e = 0; // combined thermostat: Ohmic field sample, as d_noise
         int pregen_host = 0;   // pre-generated on the host (one upload per step)

         // pre-generated on the GPU: overlap-save window state (as fft_t in quantum/internal.hpp)
         namespace fftg{
            int seg = 0, nc = 0, nk = 0, nv = 0, M = 1, n_tables = 1, lorentzian = 0;
            double norm = 0.0;
            size_t R = 0;             // realisations, 3 per atom: r = 3*atom + component
            int B = 0, B_last = 0;    // realisations per cuFFT batch, and in the last batch
            bool primed = false;
            uint64_t step = 0, window_start = 0;
            double* d_tail = 0;       // [r][2 seg] white input carried over to the next window
            double* d_coarse = 0;     // [r][4 seg] shaped coarse samples of the current window
            double* d_work = 0;       // [b][6 seg] transform input / output of one batch
            double* d_psd = 0;        // [table][3 seg + 1] filter
            cufftDoubleComplex* d_spec = 0; // [b][3 seg + 1] spectrum of one batch
            cufftHandle fwd = 0, bwd = 0, fwd_last = 0, bwd_last = 0;
            bool has_last = false;    // plans for a shorter last batch exist
         }

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
         // ohmic_bath selects the output amplitude: 0 spin (oscillator) bath, 1 Ohmic bath.
         __global__ void bath_kernel(double* s, double* out, const int* material,
                                     curandStatePhilox4_32_10_t* states,
                                     const int n_atoms, const int n_modes, const int spectrum,
                                     const double white_amp, const int ohmic_bath){
            const int tid = blockIdx.x*blockDim.x + threadIdx.x;
            curandStatePhilox4_32_10_t rng = states[tid];
            normal_source_t gauss(&rng);
            const size_t N = n_atoms;
            for(int atom = tid; atom < n_atoms; atom += blockDim.x*gridDim.x){
               const double amplitude = ohmic_bath ? c_amplitude_e[material[atom]] : c_amplitude[material[atom]];
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

         // Classical open-system thermostat: white noise into the oscillator (ohmic_bath 0)
         // or the Ohmic field of the combined thermostat (ohmic_bath 1)
         __global__ void white_kernel(double* out, const int* material,
                                      curandStatePhilox4_32_10_t* states, const int n_atoms,
                                      const int ohmic_bath){
            const int tid = blockIdx.x*blockDim.x + threadIdx.x;
            curandStatePhilox4_32_10_t rng = states[tid];
            const size_t N = n_atoms;
            for(int atom = tid; atom < n_atoms; atom += blockDim.x*gridDim.x){
               const double w = ohmic_bath ? c_white_e[material[atom]] : c_white[material[atom]];
               const float4 g = curand_normal4(&rng);
               out[atom]       = w*g.x;
               out[N + atom]   = w*g.y;
               out[2*N + atom] = w*g.z;
            }
            states[tid] = rng;
         }

         //---------------------------------------------------------------------
         // Pre-generated noise on the GPU: one window for the realisations
         // r0 ... r0 + nb - 1 (see fft_generate_window in quantum/bath_fft.cpp)
         //---------------------------------------------------------------------
         // transform input: the carried-over context (fresh noise in the first window), then fresh white noise
         __global__ void fft_fill_kernel(double* work, const double* tail, const size_t r0, const int nb,
                                         const int nc, const int overlap, const int primed,
                                         curandStatePhilox4_32_10_t* states){
            const int tid = blockIdx.x*blockDim.x + threadIdx.x;
            curandStatePhilox4_32_10_t rng = states[tid];
            normal_source_t gauss(&rng);
            const size_t total = static_cast<size_t>(nb)*nc;
            for(size_t e = tid; e < total; e += static_cast<size_t>(blockDim.x)*gridDim.x){
               const size_t b = e/nc;
               const int i = static_cast<int>(e - b*nc);
               work[e] = (primed && i < overlap) ? tail[(r0 + b)*overlap + i] : gauss.next();
            }
            states[tid] = rng;
         }

         // the last two input segments become the next window's context
         __global__ void fft_save_tail_kernel(const double* work, double* tail, const size_t r0, const int nb,
                                              const int nc, const int overlap, const int keep_from){
            const size_t total = static_cast<size_t>(nb)*overlap;
            for(size_t e = blockIdx.x*blockDim.x + threadIdx.x; e < total; e += static_cast<size_t>(blockDim.x)*gridDim.x){
               const size_t b = e/overlap;
               const int i = static_cast<int>(e - b*overlap);
               tail[(r0 + b)*overlap + i] = work[b*nc + keep_from + i];
            }
         }

         // filter: one table per material (thermostat) or one table (llg-heun bath)
         __global__ void fft_filter_kernel(cufftDoubleComplex* spec, const double* psd, const int* material,
                                           const size_t r0, const int nb, const int nk, const int lorentzian){
            const size_t total = static_cast<size_t>(nb)*nk;
            for(size_t e = blockIdx.x*blockDim.x + threadIdx.x; e < total; e += static_cast<size_t>(blockDim.x)*gridDim.x){
               const size_t b = e/nk;
               const int k = static_cast<int>(e - b*nk);
               const int table = lorentzian ? material[(r0 + b)/3] : 0;
               const double f = psd[static_cast<size_t>(table)*nk + k];
               spec[e].x *= f;
               spec[e].y *= f;
            }
         }

         // keep the middle four segments, normalised
         __global__ void fft_extract_kernel(const double* work, double* coarse, const size_t r0, const int nb,
                                            const int nc, const int nv, const int seg, const double norm){
            const size_t total = static_cast<size_t>(nb)*nv;
            for(size_t e = blockIdx.x*blockDim.x + threadIdx.x; e < total; e += static_cast<size_t>(blockDim.x)*gridDim.x){
               const size_t b = e/nv;
               const int j = static_cast<int>(e - b*nv);
               coarse[(r0 + b)*nv + j] = work[b*nc + seg + j]*norm;
            }
         }

         // this step's sample: linear interpolation between coarse samples j and j + 1, times the amplitude
         __global__ void fft_interpolate_kernel(const double* coarse, double* out, const int* material,
                                                const int n_atoms, const int nv, const int j, const double frac){
            const size_t N = n_atoms;
            for(int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < n_atoms; atom += blockDim.x*gridDim.x){
               const double a = c_amplitude[material[atom]];
               for(int c = 0; c < 3; c++){
                  const size_t base = (3*static_cast<size_t>(atom) + c)*nv;
                  out[c*N + atom] = a*(coarse[base + j]*(1.0 - frac) + coarse[base + j + 1]*frac);
               }
            }
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

         // Oscillator at its equilibrium for the current spins, q = A S/omega0^2, p = 0 (see quantum/llg.cpp)
         __global__ void seed_oscillator_kernel(const cu_real_t* sx, const cu_real_t* sy, const cu_real_t* sz,
                                                const int* material, double* q, double* p, const int n){
            for(int i = blockIdx.x*blockDim.x + threadIdx.x; i < n; i += blockDim.x*gridDim.x){
               const int m = material[i];
               const double f = c_A[m]/(c_omega0[m]*c_omega0[m]);
               q[i] = f*sx[i]; q[n + i] = f*sy[i]; q[2*n + i] = f*sz[i];
               p[i] = 0.0;     p[n + i] = 0.0;     p[2*n + i] = 0.0;
            }
         }

         // RK4: save (S, q, p)
         __global__ void save_kernel(const cu_real_t* sx, const cu_real_t* sy, const cu_real_t* sz,
                                     const double* q, const double* p, double* yin, const int n_atoms){
            const size_t N = n_atoms;
            for(int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < n_atoms; atom += blockDim.x*gridDim.x){
               yin[atom] = sx[atom]; yin[N + atom] = sy[atom]; yin[2*N + atom] = sz[atom];
               for(int c = 0; c < 3; c++){
                  yin[(3 + c)*N + atom] = q[c*N + atom];
                  yin[(6 + c)*N + atom] = p[c*N + atom];
               }
            }
         }

         // RK4 stage: slope from the source state, then the next predictor (stages 1-3),
         // whose spin is published so the spin fields can be recomputed.
         __global__ void stage_kernel(const int stage, const double dt,
                                      const double* yin, double* ypred, double* k_all,
                                      const cu_real_t* hsx, const cu_real_t* hsy, const cu_real_t* hsz,
                                      const cu_real_t* hex, const cu_real_t* hey, const cu_real_t* hez,
                                      const double* noise, const double* noise_e, const int* material,
                                      cu_real_t* sx, cu_real_t* sy, cu_real_t* sz, const int n_atoms,
                                      const int pregen){
            const size_t N = n_atoms;
            const double* ysrc = (stage == 1) ? yin : ypred;
            double* k = k_all + static_cast<size_t>(stage - 1)*9*N;
            const double coeff = (stage == 3) ? dt : 0.5*dt;
            for(int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < n_atoms; atom += blockDim.x*gridDim.x){
               double y[9], y0[9], dydt[9], H[3], xi[3];
               for(int i = 0; i < 9; i++){ y[i] = ysrc[i*N + atom]; y0[i] = yin[i*N + atom]; }
               H[0] = hsx[atom] + hex[atom];
               H[1] = hsy[atom] + hey[atom];
               H[2] = hsz[atom] + hez[atom];
               xi[0] = noise[atom]; xi[1] = noise[N + atom]; xi[2] = noise[2*N + atom];
               const int m = material[atom];
               if(pregen){
                  // pre-generated: the Lorentzian-shaped sample acts on the spin, the oscillator has no noise
                  H[0] += xi[0]; H[1] += xi[1]; H[2] += xi[2];
                  quantum::internal::fft_equation_of_motion(y, H, c_A[m], c_gamma[m], c_omega0[m], dydt);
               }
               else if(noise_e){
                  const double h_e[3] = {noise_e[atom], noise_e[N + atom], noise_e[2*N + atom]};
                  quantum::internal::ho_ohmic_equation_of_motion(y, H, xi, h_e, c_A[m], c_gamma[m], c_omega0[m], c_alpha_e[m], dydt);
               }
               else quantum::internal::ho_equation_of_motion(y, H, xi, c_A[m], c_gamma[m], c_omega0[m], dydt);
               for(int i = 0; i < 9; i++) k[i*N + atom] = dydt[i];
               if(stage < 4){
                  double yp[9];
                  quantum::internal::rk4_predict(y0, dydt, coeff, yp);
                  for(int i = 0; i < 9; i++) ypred[i*N + atom] = yp[i];
                  sx[atom] = yp[0]; sy[atom] = yp[1]; sz[atom] = yp[2];
               }
            }
         }

         // RK4: combine the four slopes and commit (S, q, p)
         __global__ void finish_kernel(const double dt_over_6, const double* yin, const double* k_all,
                                       cu_real_t* sx, cu_real_t* sy, cu_real_t* sz,
                                       double* q, double* p, const int n_atoms){
            const size_t N = n_atoms;
            for(int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < n_atoms; atom += blockDim.x*gridDim.x){
               double y0[9], k1[9], k2[9], k3[9], k4[9], y[9];
               for(int i = 0; i < 9; i++){
                  y0[i] = yin[i*N + atom];
                  k1[i] = k_all[(0*9 + i)*N + atom];
                  k2[i] = k_all[(1*9 + i)*N + atom];
                  k3[i] = k_all[(2*9 + i)*N + atom];
                  k4[i] = k_all[(3*9 + i)*N + atom];
               }
               quantum::internal::rk4_combine(y0, k1, k2, k3, k4, dt_over_6, y);
               sx[atom] = y[0]; sy[atom] = y[1]; sz[atom] = y[2];
               for(int c = 0; c < 3; c++){ q[c*N + atom] = y[3 + c]; p[c*N + atom] = y[6 + c]; }
            }
         }

         //---------------------------------------------------------------------
         // Host functions
         //---------------------------------------------------------------------
         static void cufft_check(const cufftResult r, const char* what){
            if(r != CUFFT_SUCCESS){
               err::fatal(std::string("Error - cuFFT failed (") + what + "), code " + std::to_string(static_cast<int>(r)) + ".");
            }
         }

         // Pre-generated noise on the GPU: sizes from the host setup (fft_setup), buffers, plans
         static void fft_initialise(const quantum::gpu_view_t& v){
            namespace g = fftg;
            g::seg = v.fft_seg; g::nc = 6*g::seg; g::nk = g::nc/2 + 1; g::nv = 4*g::seg; g::M = v.fft_M;
            g::n_tables = v.fft_n_tables; g::lorentzian = (v.mode == 1) ? 1 : 0; g::norm = v.fft_norm;
            g::R = 3*static_cast<size_t>(n_atoms);

            // batch of realisations per transform: about 256 MiB of work buffers
            const double per_realisation = g::nc*sizeof(double) + g::nk*sizeof(cufftDoubleComplex);
            g::B = static_cast<int>(std::min<double>(static_cast<double>(g::R), std::max(1.0, 256.0*1024.0*1024.0/per_realisation)));
            g::B_last = static_cast<int>(g::R - (g::R - 1)/g::B*g::B);

            const double window = static_cast<double>(g::R)*(g::nv + 2*g::seg)*sizeof(double);
            const double work = g::B*per_realisation + static_cast<double>(g::n_tables)*g::nk*sizeof(double);
            size_t free_bytes = 0, total_bytes = 0;
            cudaMemGetInfo(&free_bytes, &total_bytes);
            if(window + work > 0.9*free_bytes){
               std::ostringstream message;
               message << "Error - pre-generated quantum noise on the GPU needs " << (window + work)/(1024.0*1024.0*1024.0)
                       << " GiB of device memory, " << free_bytes/(1024.0*1024.0*1024.0) << " GiB are free. "
                       << "Set quantum:noise-window-size to a smaller number of steps (memory grows with it), "
                       << "or use quantum:pre-generated-device = host.";
               err::fatal(message.str());
            }

            cudaMalloc((void**)&g::d_tail, g::R*2*g::seg*sizeof(double));
            cudaMalloc((void**)&g::d_coarse, g::R*g::nv*sizeof(double));
            cudaMalloc((void**)&g::d_work, static_cast<size_t>(g::B)*g::nc*sizeof(double));
            cudaMalloc((void**)&g::d_spec, static_cast<size_t>(g::B)*g::nk*sizeof(cufftDoubleComplex));
            cudaMalloc((void**)&g::d_psd, static_cast<size_t>(g::n_tables)*g::nk*sizeof(double));
            cudaMemcpy(g::d_psd, v.fft_sqrt_psd, static_cast<size_t>(g::n_tables)*g::nk*sizeof(double), cudaMemcpyHostToDevice);

            int n[1] = {g::nc};
            cufft_check(cufftPlanMany(&g::fwd, 1, n, NULL, 1, g::nc, NULL, 1, g::nk, CUFFT_D2Z, g::B), "plan D2Z");
            cufft_check(cufftPlanMany(&g::bwd, 1, n, NULL, 1, g::nk, NULL, 1, g::nc, CUFFT_Z2D, g::B), "plan Z2D");
            g::has_last = (g::B_last != g::B);
            if(g::has_last){
               cufft_check(cufftPlanMany(&g::fwd_last, 1, n, NULL, 1, g::nc, NULL, 1, g::nk, CUFFT_D2Z, g::B_last), "plan D2Z");
               cufft_check(cufftPlanMany(&g::bwd_last, 1, n, NULL, 1, g::nk, NULL, 1, g::nc, CUFFT_Z2D, g::B_last), "plan Z2D");
            }
            check_cuda_errors(__FILE__, __LINE__);

            zlog << zTs() << "Quantum pre-generated noise on the GPU (cuFFT): " << g::R << " realisations, FFT length "
                 << g::nc << ", window " << g::nv*g::M << " steps, " << (g::R + g::B - 1)/g::B << " batches of " << g::B
                 << ", device memory " << (window + work)/(1024.0*1024.0) << " MiB" << std::endl;
         }

         // one window for every realisation, batch by batch (fft_generate_window)
         static void fft_generate_window(){
            namespace g = fftg;
            const int overlap = 2*g::seg;
            for(size_t r0 = 0; r0 < g::R; r0 += g::B){
               const int nb = static_cast<int>(std::min<size_t>(g::B, g::R - r0));
               const bool last = (nb != g::B);
               fft_fill_kernel <<< cu::grid_size, cu::block_size >>> (g::d_work, g::d_tail, r0, nb, g::nc, overlap, g::primed ? 1 : 0, d_rng);
               fft_save_tail_kernel <<< cu::grid_size, cu::block_size >>> (g::d_work, g::d_tail, r0, nb, g::nc, overlap, 4*g::seg);
               cufft_check(cufftExecD2Z(last ? g::fwd_last : g::fwd, g::d_work, g::d_spec), "D2Z");
               fft_filter_kernel <<< cu::grid_size, cu::block_size >>> (g::d_spec, g::d_psd, cu::atoms::d_materials, r0, nb, g::nk, g::lorentzian);
               cufft_check(cufftExecZ2D(last ? g::bwd_last : g::bwd, g::d_spec, g::d_work), "Z2D");
               fft_extract_kernel <<< cu::grid_size, cu::block_size >>> (g::d_work, g::d_coarse, r0, nb, g::nc, g::nv, g::seg, g::norm);
            }
            g::primed = true;
            check_cuda_errors(__FILE__, __LINE__);
         }

         // this step's sample into d_noise (fft_draw), a new window when the current one is used up
         static void fft_draw_device(){
            namespace g = fftg;
            if(!g::primed || (g::step - g::window_start) >= static_cast<uint64_t>(g::nv*g::M - g::M)){
               g::window_start = g::step;
               fft_generate_window();
            }
            const uint64_t local = g::step - g::window_start;
            const int j = static_cast<int>(local/g::M);
            const double frac = static_cast<double>(local % g::M)/g::M;
            fft_interpolate_kernel <<< cu::grid_size, cu::block_size >>> (g::d_coarse, d_noise, cu::atoms::d_materials, n_atoms, g::nv, j, frac);
            g::step++;
         }

         void initialise(){

            quantum::gpu_view_t v;
            quantum::gpu_view(v);

            n_atoms = ::atoms::num_atoms;
            n_modes = v.n_modes;
            spectrum = v.spectrum;
            mode = v.mode;
            ohmic = v.ohmic;
            pregen = v.pregenerated;
            pregen_host = v.pregenerated_host;

            if(!pregen && n_modes > max_modes) err::fatal("Error - quantum:bath-modes exceeds the GPU limit of 256.");
            if(v.n_materials > max_materials) err::fatal("Error - the quantum module on the GPU supports at most 64 materials.");

            const int n_threads = cu::grid_size*cu::block_size;
            cudaMalloc((void**)&d_rng, n_threads*sizeof(curandStatePhilox4_32_10_t));
            init_rng_kernel <<< cu::grid_size, cu::block_size >>> (d_rng, 0x9e3779b97f4a7c15ULL ^ static_cast<unsigned long long>(::mtrandom::integration_seed), n_threads);

            const size_t N = n_atoms;
            cudaMalloc((void**)&d_noise, 3*N*sizeof(double));
            cudaMemset(d_noise, 0, 3*N*sizeof(double));
            if(spectrum != 0 && !pregen){
               cudaMalloc((void**)&d_s, 3*N*n_modes*sizeof(double));
               cudaMemset(d_s, 0, 3*N*n_modes*sizeof(double));   // as the CPU bath: modes start from zero
            }
            if(mode == 1){
               cudaMalloc((void**)&d_q, 3*N*sizeof(double));     cudaMemset(d_q, 0, 3*N*sizeof(double));
               cudaMalloc((void**)&d_p, 3*N*sizeof(double));     cudaMemset(d_p, 0, 3*N*sizeof(double));
               cudaMalloc((void**)&d_yin, 9*N*sizeof(double));
               cudaMalloc((void**)&d_ypred, 9*N*sizeof(double));
               cudaMalloc((void**)&d_k, 4*9*N*sizeof(double));
               cudaMemcpyToSymbol(c_A, v.A, v.n_materials*sizeof(double));
               cudaMemcpyToSymbol(c_gamma, v.gamma, v.n_materials*sizeof(double));
               cudaMemcpyToSymbol(c_omega0, v.omega0, v.n_materials*sizeof(double));
               if(ohmic){
                  cudaMalloc((void**)&d_noise_e, 3*N*sizeof(double));
                  cudaMemset(d_noise_e, 0, 3*N*sizeof(double));
                  if(spectrum != 0){
                     cudaMalloc((void**)&d_s_e, 3*N*n_modes*sizeof(double));
                     cudaMemset(d_s_e, 0, 3*N*n_modes*sizeof(double));
                  }
                  cudaMemcpyToSymbol(c_alpha_e, v.alpha_e, v.n_materials*sizeof(double));
                  cudaMemcpyToSymbol(c_amplitude_e, v.amplitude_e, v.n_materials*sizeof(double)); // same in both phases
               }
               // first step: the device spins are final (also after a checkpoint load)
               seed_oscillator_kernel <<< cu::grid_size, cu::block_size >>> (cu::atoms::d_x_spin, cu::atoms::d_y_spin,
                  cu::atoms::d_z_spin, cu::atoms::d_materials, d_q, d_p, n_atoms);
            }
            if(pregen && !pregen_host) fft_initialise(v);
            check_cuda_errors(__FILE__, __LINE__);

            zlog << zTs() << "Quantum module on the GPU: " << (mode == 1 ? (ohmic ? "open-system LLG + Ohmic bath (RK4)" : "open-system LLG (RK4)") : "coloured bath for llg-heun")
                 << ", spectrum " << (spectrum == 0 ? "classical" : spectrum == 1 ? "quantum" : "quantum-no-zero")
                 << (pregen ? (pregen_host ? ", pre-generated (FFT on the host, one upload per step)" : ", pre-generated (cuFFT on the GPU)")
                            : (spectrum == 0 ? ", white noise" : spectrum == 1 ? ", Matsubara OU modes" : ", log-spaced OU modes"))
                 << ", " << (pregen ? 0 : n_modes) << " modes, " << n_atoms << " atoms" << std::endl;

            initialised = true;
            return;
         }

         // Upload coefficients when the temperature or the phase changed; return the view
         void prepare(quantum::gpu_view_t& v){
            if(!initialised) initialise();
            quantum::gpu_view(v);
            if(!pregen && v.spectrum != 0 && v.T_scaled != uploaded_T){
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
            if(mode == 1 && v.spectrum == 0 && !pregen){
               // the classical prefactor carries the temperature: refresh every step (n_materials doubles)
               cudaMemcpyToSymbol(c_white, v.white.data(), v.n_materials*sizeof(double));
               if(ohmic) cudaMemcpyToSymbol(c_white_e, v.white_e.data(), v.n_materials*sizeof(double));
            }
            return;
         }

         void draw(const quantum::gpu_view_t& v){
            // pre-generated: the CPU generator draws on the host, one upload per step
            // ([component][atom], as the device baths write it)
            if(pregen && !pregen_host){
               fft_draw_device();
               check_cuda_errors(__FILE__, __LINE__);
               return;
            }
            if(pregen){
               const double* h_x = 0;
               const double* h_y = 0;
               const double* h_z = 0;
               quantum::gpu_pregenerated_draw(h_x, h_y, h_z);
               const size_t N = n_atoms;
               cudaMemcpy(d_noise,         h_x, N*sizeof(double), cudaMemcpyHostToDevice);
               cudaMemcpy(d_noise + N,     h_y, N*sizeof(double), cudaMemcpyHostToDevice);
               cudaMemcpy(d_noise + 2*N,   h_z, N*sizeof(double), cudaMemcpyHostToDevice);
               check_cuda_errors(__FILE__, __LINE__);
               return;
            }
            // the Ohmic bath of the combined thermostat shares the coefficients, not the samples
            for(int e = 0; e <= ohmic; e++){
               double* out = e ? d_noise_e : d_noise;
               if(v.spectrum == 0){
                  white_kernel <<< cu::grid_size, cu::block_size >>> (out, cu::atoms::d_materials, d_rng, n_atoms, e);
               }
               else{
                  bath_kernel <<< cu::grid_size, cu::block_size >>> (e ? d_s_e : d_s, out, cu::atoms::d_materials, d_rng,
                                                                      n_atoms, n_modes, v.spectrum, v.white_amp, e);
               }
            }
            check_cuda_errors(__FILE__, __LINE__);
         }

         // quantum:export-noise: copy the field of the export atom ([component][atom]
         // array on the device) to the host and append it to the export file;
         // add_sample: plus this step's pre-generated sample in d_noise (thermostat:
         // q + sample, as export_thermostat_sample() on the CPU)
         void export_field(const double* d_field, const bool add_sample = false){
            const int a = quantum::gpu_export_atom();
            if(a < 0 || a >= n_atoms) return;
            const size_t N = n_atoms;
            double v[3];
            for(int c = 0; c < 3; c++){
               cudaMemcpy(&v[c], d_field + c*N + a, sizeof(double), cudaMemcpyDeviceToHost);
               if(add_sample){
                  double sample = 0.0;
                  cudaMemcpy(&sample, d_noise + c*N + a, sizeof(double), cudaMemcpyDeviceToHost);
                  v[c] += sample;
               }
            }
            quantum::gpu_export_sample(v[0], v[1], v[2]);
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
         qg::export_field(qg::d_noise);    // as quantum::generate() on the CPU: the bath sample of this step
      }

      void quantum_finalize(){
         namespace qg = quantum_gpu;
         if(!qg::initialised) return;
         cudaFree(qg::d_rng); cudaFree(qg::d_noise);
         if(qg::d_s) cudaFree(qg::d_s);
         if(qg::d_q){ cudaFree(qg::d_q); cudaFree(qg::d_p); cudaFree(qg::d_yin); cudaFree(qg::d_ypred); cudaFree(qg::d_k); }
         if(qg::d_noise_e) cudaFree(qg::d_noise_e);
         if(qg::d_s_e) cudaFree(qg::d_s_e);
         if(qg::fftg::d_coarse){
            namespace g = qg::fftg;
            cudaFree(g::d_tail); cudaFree(g::d_coarse); cudaFree(g::d_work); cudaFree(g::d_spec); cudaFree(g::d_psd);
            cufftDestroy(g::fwd); cufftDestroy(g::bwd);
            if(g::has_last){ cufftDestroy(g::fwd_last); cufftDestroy(g::bwd_last); }
            g::d_coarse = 0;
         }
         qg::initialised = false;
      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Open-system LLG step on the GPU (mirrors quantum::internal::llg_ho_serial)
   //---------------------------------------------------------------------------
   void quantum_llg(){

      namespace qg = cu::quantum_gpu;

      quantum::gpu_view_t v;
      qg::prepare(v);
      const int n = qg::n_atoms;
      const double dt = ::mp::dt;

      qg::save_kernel <<< cu::grid_size, cu::block_size >>> (cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
                                                             qg::d_q, qg::d_p, qg::d_yin, n);
      cu::update_spin_fields();
      cu::update_external_fields();     // applied (+FMR, dipole) fields only: the thermal part is off in this mode
      qg::draw(v);                      // one bath sample per step, seen by all four stages (pre-generated: from the host)

      for(int stage = 1; stage <= 4; stage++){
         qg::stage_kernel <<< cu::grid_size, cu::block_size >>> (stage, dt, qg::d_yin, qg::d_ypred, qg::d_k,
               cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
               cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
               qg::d_noise, qg::ohmic ? qg::d_noise_e : nullptr, cu::atoms::d_materials,
               cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin, n, qg::pregen);
         if(stage < 4) cu::update_spin_fields();
      }

      qg::finish_kernel <<< cu::grid_size, cu::block_size >>> (dt/6.0, qg::d_yin, qg::d_k,
            cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin, qg::d_q, qg::d_p, n);
      check_cuda_errors(__FILE__, __LINE__);
      qg::export_field(qg::d_q, qg::pregen); // as export_thermostat_sample() on the CPU: q (+ sample) after the step

      return;
   }

} // end of vcuda namespace

#endif
