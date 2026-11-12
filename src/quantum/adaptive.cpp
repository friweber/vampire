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
//   Adaptive quantum bath (quantum:adaptive-bath), the adQTB of molecular
//   dynamics (Mangaud et al. 2019, Mauger et al. 2021) for the LLG equation.
//
//   The pre-generated noise spectrum gets a gain g(w) >= 0 per frequency band,
//   shared by all materials: C_bb(w) = g(w) C_target(w). Before every new noise
//   window the probe series of the FDT monitor (fdt_monitor.cpp; block means of
//   b and ds/dt since the last update) are transformed and
//      P_in(w)   = mu Re C_{b sdot}(w)          (power the noise injects)
//      P_diss(w) = (alpha mu/gamma) C_{sdot sdot}(w)  (power the damping removes)
//   are summed over probes and components and binned into bands. The band holds
//   rho = g P_diss/P_in of its target energy (FDT: rho = 1); the gain is moved
//   multiplicatively, g <- g rho^-eta, by at most 20 % per update, using spectra
//   averaged over the updates with weight 1/2 per update (exponential memory).
//   Bands without measurable injected power keep their gain. The tables of the
//   noise filter become sqrt_psd0 sqrt(g(w)), g linearly interpolated between
//   the band centres; the GPU re-uploads them (table_version).
//
//   Output: adaptive_gain.dat (one line per update: time, then g and rho per
//   band; header with the band centres in meV) and the final gain table
//   quantum:adaptive-gain-file-out (default adaptive_gain_final.dat), which a
//   later run can start from (quantum:adaptive-gain-file-in).
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

// Vampire headers
#include "constants.hpp"
#include "errors.hpp"
#include "integrate.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vio.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      bool adaptive_bath = false;
      double adaptive_eta = 0.1;
      int adaptive_bands = 32;
      double adaptive_max_energy_meV = 0.0;
      int64_t adaptive_freeze_step = -1;
      std::string adaptive_gain_in = "";
      std::string adaptive_gain_out = "adaptive_gain_final.dat";

      namespace adapt{
         bool ready = false;
         std::vector<double> w_centre;   // band centres (rad/s)
         double w_top = 0.0;             // upper edge of the last band (rad/s)
         std::vector<double> gain;       // g per band
         std::vector<double> acc_in, acc_d, rho;
         std::ofstream log;
         int n_updates = 0;

         // band layout: uniform on [0, w_top]; w_top = noise cut-off, else the Nyquist frequency of the block series
         void layout(){
            const double hbar = constants::hbar, meV = 1.0e-3*constants::eV;
            const double w_nyq = M_PI/fdt::block_dt_SI;
            if(adaptive_max_energy_meV > 0.0) w_top = std::min(w_nyq, adaptive_max_energy_meV*meV/hbar);
            else if(noise_cutoff_energy_meV > 0.0) w_top = std::min(w_nyq, 1.2*noise_cutoff_energy_meV*meV/hbar);
            else err::fatal("Error - quantum:adaptive-bath needs quantum:adaptive-max-energy (meV) or quantum:noise-cutoff-energy "
                            "to place its frequency bands.");
            const int n = adaptive_bands;
            w_centre.resize(n);
            for(int i = 0; i < n; i++) w_centre[i] = (i + 0.5)*w_top/n;
            if(gain.size() != static_cast<size_t>(n)) gain.assign(n, 1.0);
            acc_in.assign(n, 0.0); acc_d.assign(n, 0.0); rho.assign(n, std::nan(""));
            log.open("adaptive_gain.dat");
            log << "# VAMPIRE adaptive bath: gain g and band energy rho = g P_diss/P_in per band; eta " << adaptive_eta
                << ", block " << fdt::block_dt_SI*1e15 << " fs, " << fdt::probe.size() << " probe atoms\n# band centres (meV):";
            for(int i = 0; i < n; i++) log << " " << w_centre[i]*hbar/meV;
            log << "\n# time(ps) n_blocks then g_1..g_n then rho_1..rho_n\n";
            ready = true;
         }

         // real FFT of a Hann-windowed, mean-free series; bins 0..n/2 (FFTW, plan cached per length)
         void spectrum(const std::vector<double>& x, std::vector<double>& re, std::vector<double>& im){
            const int n = x.size(); const int nf = n/2;
            re.assign(nf + 1, 0.0); im.assign(nf + 1, 0.0);
            double mean = 0.0; for(double v : x) mean += v; mean /= n;
            #ifdef FFT
            static int plan_n = 0; static fftw_plan plan = 0; static double* in = 0; static fftw_complex* out = 0;
            if(plan_n != n){
               if(plan){ fftw_destroy_plan(plan); fftw_free(in); fftw_free(out); }
               in = static_cast<double*>(fftw_malloc(sizeof(double)*n));
               out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*(nf + 1)));
               plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
               plan_n = n;
            }
            for(int t = 0; t < n; t++) in[t] = (x[t] - mean)*(0.5 - 0.5*std::cos(2.0*M_PI*t/(n - 1)));
            fftw_execute(plan);
            for(int k = 0; k <= nf; k++){ re[k] = out[k][0]; im[k] = out[k][1]; }
            #else
            err::fatal("Error - quantum:adaptive-bath requires FFTW (-DFFT -lfftw3).");
            #endif
         }

      } // end of adapt namespace

      //---------------------------------------------------------------------------
      // gain from file (two columns: band centre in meV, g), interpolated onto
      // the bands of this run; applied to the tables right away
      //---------------------------------------------------------------------------
      void adaptive_setup(bath_t& b){
         if(!adaptive_bath) return;
         if(integrate::integrator != integrate::integrator_t::llg_heun){
            err::fatal("Error - quantum:adaptive-bath works with sim:integrator = llg-heun only.");
         }
         fdt::ensure_initialised();
         adapt::layout();
         if(!adaptive_gain_in.empty()){
            std::ifstream f(adaptive_gain_in.c_str());
            if(!f.is_open()) err::fatal("Error - cannot read quantum:adaptive-gain-file-in " + adaptive_gain_in);
            std::vector<double> e, g; std::string line;
            while(std::getline(f, line)){
               if(line.empty() || line[0] == '#') continue;
               std::istringstream is(line); double a, c;
               if(is >> a >> c){ e.push_back(a); g.push_back(c); }
            }
            const double meV = 1.0e-3*constants::eV;
            for(size_t i = 0; i < adapt::gain.size(); i++){
               const double x = adapt::w_centre[i]*constants::hbar/meV;
               if(e.empty()) break;
               if(x <= e.front()) adapt::gain[i] = g.front();
               else if(x >= e.back()) adapt::gain[i] = g.back();
               else{
                  size_t j = std::upper_bound(e.begin(), e.end(), x) - e.begin();
                  adapt::gain[i] = g[j-1] + (g[j] - g[j-1])*(x - e[j-1])/(e[j] - e[j-1]);
               }
            }
            zlog << zTs() << "Adaptive bath: gain read from " << adaptive_gain_in << std::endl;
         }
         fft_apply_gain(b, adapt::w_centre, adapt::gain);
         zlog << zTs() << "Adaptive bath: " << adaptive_bands << " bands up to " << adapt::w_top*constants::hbar/(1e-3*constants::eV)
              << " meV, eta " << adaptive_eta << std::endl;
      }

      //---------------------------------------------------------------------------
      // update before a new window; true if the tables changed
      //---------------------------------------------------------------------------
      bool adaptive_before_window(bath_t& b){
         if(!adaptive_bath || !adapt::ready) return false;
         if(adaptive_freeze_step >= 0 && static_cast<int64_t>(sim::time) >= adaptive_freeze_step){
            fdt::series_b.clear(); fdt::series_v.clear(); fdt::series_blocks = 0;
            return false;
         }
         const int nblk = fdt::series_blocks;
         // keep collecting (over several windows if needed) until every band holds >= 3 spectral bins:
         // 2 pi/(nblk dt_block) <= band width/3
         const double band_width = adapt::w_top/adaptive_bands;
         const int nblk_min = std::max(32, static_cast<int>(std::ceil(3.0*2.0*M_PI/(fdt::block_dt_SI*band_width))));
         if(nblk < nblk_min) return false;
         const int np = fdt::probe.size();
         const int nf = nblk/2;
         const double dw = 2.0*M_PI/(nblk*fdt::block_dt_SI);
         const int n = adaptive_bands;
         std::vector<double> pin(n, 0.0), pd(n, 0.0);
         std::vector<double> x(nblk), y(nblk), br, bi, vr, vi;
         for(int p = 0; p < np; p++){
            for(int c = 0; c < 3; c++){
               for(int t = 0; t < nblk; t++){
                  x[t] = fdt::series_b[(static_cast<size_t>(t)*np + p)*3 + c];
                  y[t] = fdt::series_v[(static_cast<size_t>(t)*np + p)*3 + c];
               }
               adapt::spectrum(x, br, bi); adapt::spectrum(y, vr, vi);
               for(int k = 1; k <= nf; k++){
                  const double w = k*dw;
                  if(w >= adapt::w_top) break;
                  const int band = std::min(n - 1, static_cast<int>(w/adapt::w_top*n));
                  pin[band] += fdt::probe_mu[p]*(br[k]*vr[k] + bi[k]*vi[k]);
                  pd[band]  += fdt::probe_mu[p]*fdt::probe_alpha_over_gamma[p]*(vr[k]*vr[k] + vi[k]*vi[k]);
               }
            }
         }
         fdt::series_b.clear(); fdt::series_v.clear(); fdt::series_blocks = 0;

         double pmax = 0.0;
         for(int i = 0; i < n; i++){
            adapt::acc_in[i] = 0.5*adapt::acc_in[i] + pin[i];
            adapt::acc_d[i]  = 0.5*adapt::acc_d[i] + pd[i];
            pmax = std::max(pmax, adapt::acc_in[i]);
         }
         for(int i = 0; i < n; i++){
            if(adapt::acc_in[i] > 1.0e-3*pmax && adapt::acc_d[i] > 0.0){
               adapt::rho[i] = adapt::gain[i]*adapt::acc_d[i]/adapt::acc_in[i];
               const double step = std::min(1.2, std::max(1.0/1.2, std::pow(adapt::rho[i], -adaptive_eta)));
               adapt::gain[i] = std::min(100.0, std::max(1.0e-3, adapt::gain[i]*step));
            }
            else adapt::rho[i] = std::nan("");
         }
         adapt::n_updates++;
         adapt::log << sim::time*mp::dt_SI*1e12 << " " << nblk;
         for(int i = 0; i < n; i++) adapt::log << " " << adapt::gain[i];
         for(int i = 0; i < n; i++) adapt::log << " " << adapt::rho[i];
         adapt::log << "\n"; adapt::log.flush();

         fft_apply_gain(b, adapt::w_centre, adapt::gain);
         return true;
      }

      void adaptive_finalize(){
         if(!adaptive_bath || !adapt::ready) return;
         std::ofstream f(adaptive_gain_out.c_str());
         const double meV = 1.0e-3*constants::eV;
         f << "# VAMPIRE adaptive bath, final gain after " << adapt::n_updates << " updates: band centre (meV), g, last rho\n";
         for(size_t i = 0; i < adapt::gain.size(); i++)
            f << adapt::w_centre[i]*constants::hbar/meV << " " << adapt::gain[i] << " " << adapt::rho[i] << "\n";
         adapt::log.close();
         zlog << zTs() << "Adaptive bath: " << adapt::n_updates << " updates, final gain in " << adaptive_gain_out << std::endl;
      }

   } // end of internal namespace

} // end of quantum namespace
