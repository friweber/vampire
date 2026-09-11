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
#include <limits>

// Vampire headers
#include "quantum.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// quantum module headers
#include "internal.hpp"
#include "llg_atom.hpp"

namespace quantum{

   namespace internal{

      //---------------------------------------------------------------------------
      // Squared error between the bath spectrum of n_modes Ornstein-Uhlenbeck
      // modes log-spaced in [10^log_lo, 10^log_hi] and the target
      // w (coth(w/2T) - 1), summed over the frequency grid
      //---------------------------------------------------------------------------
      static double fit_error(const double log_lo, const double log_hi, const double T_scaled,
                              const int n_modes, const std::vector<double>& omega){

         if(log_hi <= log_lo + 0.05) return std::numeric_limits<double>::infinity();

         std::vector<double> lambda(n_modes);
         for(int k = 0; k < n_modes; k++){
            lambda[k] = std::pow(10.0, log_lo + static_cast<double>(k)/(n_modes - 1)*(log_hi - log_lo));
         }

         double total = 0.0;
         for(size_t i = 0; i < omega.size(); i++){
            const double w = omega[i];
            double model = 0.0;
            for(int k = 0; k < n_modes - 1; k++){
               const double dl = lambda[k + 1] - lambda[k];
               const double x = lambda[k]/(2.0*T_scaled);
               const double coth = (std::fabs(x) < 1e-10) ? 1.0/x : 1.0/std::tanh(x);
               model += 2.0*lambda[k]*lambda[k]*(coth - 1.0)*dl/(w*w + lambda[k]*lambda[k]);
            }
            const double xw = w/(2.0*T_scaled);
            const double cothw = (std::fabs(xw) < 1e-10) ? 1.0/xw : 1.0/std::tanh(xw);
            const double target = w*(cothw - 1.0);
            total += (model - target)*(model - target);
         }

         return total;

      }

      //---------------------------------------------------------------------------
      // Function to choose the rate window of the log-spaced bath by a grid scan
      // of (log lambda_lo, log lambda_hi). The target spectrum has the
      // temperature as its only scale, so the search box is centred on T and
      // the frequency grid runs up to 15 T. The window is kept fixed when the
      // temperature moves later; only the amplitudes are refreshed.
      //---------------------------------------------------------------------------
      void log_bath_fit(bath_t& b, const double T_scaled){

         const int n_modes = b.n_modes;

         std::vector<double> omega(bath_scan_omega_points);
         const double omega_max = 15.0*T_scaled;
         for(int i = 0; i < bath_scan_omega_points; i++){
            omega[i] = (i + 1)*omega_max/bath_scan_omega_points;
         }

         const double decades = bath_scan_decades;
         const double lo_min = std::log10(T_scaled) - decades;
         const double lo_max = std::log10(T_scaled) + 1.0;
         const double hi_min = std::log10(T_scaled) - 0.5;
         const double hi_max = std::log10(T_scaled) + decades;

         const int n_scan = bath_scan_resolution;
         double best_error = std::numeric_limits<double>::infinity();
         double best_lo = lo_min;
         double best_hi = hi_min;

         for(int i = 0; i < n_scan; i++){
            const double log_lo = lo_min + i*(lo_max - lo_min)/(n_scan - 1);
            for(int j = 0; j < n_scan; j++){
               const double log_hi = hi_min + j*(hi_max - hi_min)/(n_scan - 1);
               if(log_hi <= log_lo + 0.05) continue;
               const double error = fit_error(log_lo, log_hi, T_scaled, n_modes, omega);
               if(error < best_error){
                  best_error = error;
                  best_lo = log_lo;
                  best_hi = log_hi;
               }
            }
         }

         b.lb_lambda.resize(n_modes);
         for(int k = 0; k < n_modes; k++){
            b.lb_lambda[k] = std::pow(10.0, best_lo + static_cast<double>(k)/(n_modes - 1)*(best_hi - best_lo));
         }

         if(vmpi::my_rank == 0){
            zlog << zTs() << "Quantum log-bath fit: " << n_modes << " modes, lambda in ["
                 << std::pow(10.0, best_lo) << ", " << std::pow(10.0, best_hi)
                 << "], squared error " << best_error << std::endl;
         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to build the log-bath coefficients for a given temperature.
      // Each mode carries the weight of its slice of the target spectrum,
      // sqrt(lambda_k (coth(lambda_k/2T) - 1) dlambda_k), renormalised so the
      // bath reproduces the target at low frequency. decay and coeff are the
      // exact Ornstein-Uhlenbeck update over one step.
      //---------------------------------------------------------------------------
      void log_bath_coefficients(bath_t& b, const double T_scaled){

         const int n_modes = b.n_modes;
         const std::vector<double>& lambda = b.lb_lambda;

         b.lb_decay.resize(n_modes);
         b.lb_coeff.resize(n_modes);
         b.lb_amp.resize(n_modes);

         std::vector<double> raw(n_modes, 0.0);
         for(int k = 0; k < n_modes; k++){
            const double dl = (k < n_modes - 1) ? lambda[k + 1] - lambda[k] : lambda[k] - lambda[k - 1];
            const double x = lambda[k]/(2.0*T_scaled);
            const double coth = (std::fabs(x) < 1e-10) ? 1.0/x : 1.0/std::tanh(x);
            const double weight = coth - 1.0;
            raw[k] = (weight > 0.0) ? std::sqrt(lambda[k]*weight*dl) : 0.0;
         }

         // renormalise at a reference frequency well below T
         const double omega_ref = 1.0e-3*T_scaled;
         double model = 0.0;
         for(int k = 0; k < n_modes; k++){
            model += raw[k]*raw[k]*2.0*lambda[k]/(omega_ref*omega_ref + lambda[k]*lambda[k]);
         }
         const double x_ref = omega_ref/(2.0*T_scaled);
         const double target = omega_ref*(1.0/std::tanh(x_ref) - 1.0);
         const double renorm = (model > 0.0) ? std::sqrt(target/model) : 1.0;

         for(int k = 0; k < n_modes; k++){
            const double decay = std::exp(-lambda[k]*b.dt);
            b.lb_decay[k] = decay;
            b.lb_coeff[k] = std::sqrt(1.0 - decay*decay);
            b.lb_amp[k] = raw[k]*renorm;
         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to draw one site's sample from the log-spaced bath
      //---------------------------------------------------------------------------
      void log_bath_draw(bath_t& b, const int site, const double amplitude){

         const int n_modes = b.n_modes;
         const size_t base = static_cast<size_t>(site)*n_modes;

         for(int i = 0; i < 3*n_modes; i++) b.scratch[i] = mtrandom::gaussian();

         double out[3];
         log_bath_update_atom(&b.s_x[base], &b.s_y[base], &b.s_z[base], n_modes,
                              &b.lb_decay[0], &b.lb_coeff[0], &b.lb_amp[0],
                              &b.scratch[0], out);

         b.x[site] = amplitude*out[0];
         b.y[site] = amplitude*out[1];
         b.z[site] = amplitude*out[2];

         return;

      }

   } // end of internal namespace

} // end of quantum namespace
