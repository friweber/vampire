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
#include <cmath>
#include <iostream>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      // overlap-save bookkeeping: six segments per transform, four kept
      const int segments = 6;
      const int valid_segments = 4;

      //---------------------------------------------------------------------------
      // Spectral density of the bath at frequency w and reduced temperature T.
      // All three tend to the classical plateau 2T as w -> 0; the no-zero-point
      // form is written as 2w/expm1(w/T) to keep its precision where coth - 1
      // is tiny.
      //---------------------------------------------------------------------------
      static double spectral_density(const spectrum_t spectrum, const double w, const double T){

         if(spectrum == classical || w <= 0.0) return 2.0*T;

         const double x = w/(2.0*T);

         if(spectrum == quantum_zero){
            const double coth = (x < 1e-10) ? 1.0/x : 1.0/std::tanh(x);
            return w*coth;
         }

         return (x > 350.0) ? 0.0 : 2.0*w/std::expm1(2.0*x);

      }

      //---------------------------------------------------------------------------
      // Response of the thermostat's auxiliary oscillator for material m,
      // A Gamma/((omega0^2 - w^2)^2 + Gamma^2 w^2)
      //---------------------------------------------------------------------------
      static double lorentzian_filter(const int m, const double w){

         const double omega0 = material_omega0[m];
         const double Gamma = material_gamma[m];
         const double d = (omega0*omega0 - w*w)*(omega0*omega0 - w*w) + Gamma*Gamma*w*w;
         return material_A[m]*Gamma/((d > 1e-12) ? d : 1e-12);

      }

      //---------------------------------------------------------------------------
      // Smallest integer >= n whose prime factors are all 2, 3, 5 or 7. FFTW
      // transforms of such lengths run in O(n log n); a large prime factor,
      // which the whole-run window would otherwise produce for most run
      // lengths, is far slower to plan and to execute.
      //---------------------------------------------------------------------------
      static int smooth_length(const int n){

         for(int candidate = n; ; candidate++){
            int r = candidate;
            while(r % 2 == 0) r /= 2;
            while(r % 3 == 0) r /= 3;
            while(r % 5 == 0) r /= 5;
            while(r % 7 == 0) r /= 7;
            if(r == 1) return candidate;
         }

      }

      //---------------------------------------------------------------------------
      // Function to size the window, allocate the buffers and tabulate the
      // filter sqrt(S(w) L(w)) on the window's frequency grid. With no window
      // size given the whole run is covered by one window, which needs an FFT
      // length of 6/4 of the run because two of six segments are discarded.
      // The temperature is frozen at its value here.
      //---------------------------------------------------------------------------
      void fft_setup(bath_t& b, const uint64_t n_run){

         #ifdef FFT

         fft_t& f = b.fft;
         f.M = interpolation_factor;

         // window in coarse samples
         if(window_size == 0){
            // one window for the run: interpolation reads one sample ahead and
            // the regeneration rule keeps M steps spare
            const uint64_t n_valid_needed = n_run/f.M + 3;
            f.seg = static_cast<int>((n_valid_needed + valid_segments - 1)/valid_segments);
         }
         else{
            int n_coarse_requested = static_cast<int>((window_size - 1)/f.M + 1);
            if(n_coarse_requested % segments != 0) n_coarse_requested = (n_coarse_requested/segments + 1)*segments;
            f.seg = n_coarse_requested/segments;
         }
         f.seg = smooth_length(f.seg);   // rounded up for a fast transform
         f.n_coarse = segments*f.seg;
         f.n_valid = valid_segments*f.seg;

         // refuse silently unbounded memory; the user can shrink the window
         const size_t n_realisations = 3*static_cast<size_t>(b.n_sites);
         const double bytes = static_cast<double>(n_realisations)*(f.n_valid + 2*f.seg)*sizeof(double);
         if(bytes > 8.0*1024.0*1024.0*1024.0){
            terminaltextcolor(RED);
            std::cerr << "Error: pre-generated quantum noise for " << b.n_sites << " atoms and "
                      << n_run << " steps needs " << bytes/(1024.0*1024.0*1024.0)
                      << " GiB. Set quantum:noise-window-size to a smaller number of steps, or "
                      << "use quantum:noise-generator = on-the-fly." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: pre-generated quantum noise would need " << bytes/(1024.0*1024.0*1024.0) << " GiB" << std::endl;
            err::vexit();
         }

         f.tail.assign(n_realisations*2*f.seg, 0.0);
         f.coarse.assign(n_realisations*f.n_valid, 0.0);

         f.in = static_cast<double*>(fftw_malloc(sizeof(double)*f.n_coarse));
         f.out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*(f.n_coarse/2 + 1)));
         f.result = static_cast<double*>(fftw_malloc(sizeof(double)*f.n_coarse));
         f.forward = fftw_plan_dft_r2c_1d(f.n_coarse, f.in, f.out, FFTW_MEASURE);
         f.backward = fftw_plan_dft_c2r_1d(f.n_coarse, f.out, f.result, FFTW_MEASURE);

         // filter tables: one per material when the Lorentzian is included
         const double T_scaled = sim::temperature*b.T_scale;
         const double dt_coarse = f.M*b.dt;
         const int n_tables = b.lorentzian ? static_cast<int>(material_A.size()) : 1;
         f.sqrt_psd.assign(n_tables, std::vector<double>(f.n_coarse/2 + 1, 0.0));
         for(int m = 0; m < n_tables; m++){
            for(int k = 0; k <= f.n_coarse/2; k++){
               const double w = 2.0*M_PI*k/(f.n_coarse*dt_coarse);
               double psd = spectral_density(b.spectrum, w, T_scaled);
               if(b.lorentzian) psd *= lorentzian_filter(m, w);
               f.sqrt_psd[m][k] = std::sqrt(psd);
            }
         }

         // white input of unit variance per fine step, transform normalisation,
         // and the coarse-grid variance correction, applied together on output
         f.norm = 1.0/(f.n_coarse*std::sqrt(f.M*b.dt));

         f.primed = false;
         f.step = 0;
         f.window_start = 0;

         if(vmpi::my_rank == 0){
            zlog << zTs() << "Quantum pre-generated noise: window of " << f.n_valid << " coarse samples ("
                 << f.n_valid*f.M << " steps, FFT length " << f.n_coarse << "), buffer "
                 << bytes/(1024.0*1024.0) << " MiB" << std::endl;
         }

         #else
         terminaltextcolor(RED);
         std::cerr << "Error: pre-generated quantum noise requires FFTW. Recompile with -DFFT -lfftw3 "
                   << "or use quantum:noise-generator = on-the-fly." << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
         #endif

         return;

      }

      //---------------------------------------------------------------------------
      // Function to generate one window for every realisation. The transform
      // input is the two carried-over segments followed by four segments of
      // fresh white noise; after shaping, the middle four segments are kept
      // and the last two input segments become the next window's context.
      //---------------------------------------------------------------------------
      void fft_generate_window(bath_t& b){

         #ifdef FFT

         fft_t& f = b.fft;
         const int nc = f.n_coarse;
         const int seg = f.seg;
         const int overlap = 2*seg;
         const size_t n_realisations = 3*static_cast<size_t>(b.n_sites);

         for(size_t r = 0; r < n_realisations; r++){

            const size_t tail_base = r*overlap;

            // the first window has no previous context and starts from fresh noise
            if(!f.primed){
               for(int i = 0; i < overlap; i++) f.tail[tail_base + i] = mtrandom::gaussian();
            }

            for(int i = 0; i < overlap; i++) f.in[i] = f.tail[tail_base + i];
            for(int i = overlap; i < nc; i++) f.in[i] = mtrandom::gaussian();
            for(int i = 0; i < overlap; i++) f.tail[tail_base + i] = f.in[valid_segments*seg + i];

            const int table = b.lorentzian ? atoms::type_array[r/3] : 0;
            const std::vector<double>& sqrt_psd = f.sqrt_psd[table];

            fftw_execute(f.forward);
            for(int k = 0; k <= nc/2; k++){
               f.out[k][0] *= sqrt_psd[k];
               f.out[k][1] *= sqrt_psd[k];
            }
            fftw_execute(f.backward);

            const size_t coarse_base = r*f.n_valid;
            for(int j = 0; j < f.n_valid; j++){
               f.coarse[coarse_base + j] = f.result[seg + j]*f.norm;
            }

         }

         f.primed = true;

         #endif

         return;

      }

      //---------------------------------------------------------------------------
      // Function to sample this step from the current window, generating a
      // new window on the first call and whenever the current one is used up.
      // The first window is therefore drawn after the generator is seeded.
      //---------------------------------------------------------------------------
      void fft_draw(bath_t& b, const std::vector<double>& amplitude){

         fft_t& f = b.fft;

         if(!f.primed || (f.step - f.window_start) >= static_cast<uint64_t>(f.n_valid*f.M - f.M)){
            f.window_start = f.step;
            fft_generate_window(b);
         }

         const uint64_t local = f.step - f.window_start;
         const size_t j = static_cast<size_t>(local/f.M);
         const double frac = static_cast<double>(local % f.M)/f.M;

         for(int site = 0; site < b.n_sites; site++){
            const double a = amplitude[atoms::type_array[site]];
            const size_t bx = (3*static_cast<size_t>(site) + 0)*f.n_valid;
            const size_t by = (3*static_cast<size_t>(site) + 1)*f.n_valid;
            const size_t bz = (3*static_cast<size_t>(site) + 2)*f.n_valid;
            b.x[site] = a*(f.coarse[bx + j]*(1.0 - frac) + f.coarse[bx + j + 1]*frac);
            b.y[site] = a*(f.coarse[by + j]*(1.0 - frac) + f.coarse[by + j + 1]*frac);
            b.z[site] = a*(f.coarse[bz + j]*(1.0 - frac) + f.coarse[bz + j + 1]*frac);
         }

         f.step++;

         return;

      }

      //---------------------------------------------------------------------------
      // Function to free the FFTW plans and buffers of a bath
      //---------------------------------------------------------------------------
      void fft_release(bath_t& b){

         #ifdef FFT
         fft_t& f = b.fft;
         if(f.forward){ fftw_destroy_plan(f.forward); f.forward = 0; }
         if(f.backward){ fftw_destroy_plan(f.backward); f.backward = 0; }
         if(f.in){ fftw_free(f.in); f.in = 0; }
         if(f.out){ fftw_free(f.out); f.out = 0; }
         if(f.result){ fftw_free(f.result); f.result = 0; }
         f.sqrt_psd.clear();
         f.tail.clear();
         f.coarse.clear();
         #endif

         return;

      }

   } // end of internal namespace

} // end of quantum namespace
