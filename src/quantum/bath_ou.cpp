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

// Vampire headers
#include "quantum.hpp"
#include "random.hpp"

// quantum module headers
#include "internal.hpp"
#include "llg_atom.hpp"

namespace quantum{

   namespace internal{

      //---------------------------------------------------------------------------
      // Function to build the Matsubara Ornstein-Uhlenbeck coefficients for the
      // quantum_zero spectrum w coth(w/2T) at a given reduced temperature.
      //
      // Mode n relaxes at the Matsubara frequency nu_n = 2 pi n T. Over a step
      // dt it is updated exactly, s_n <- decay_n s_n + diffuse_n xi_n, and the
      // noise it delivers is the integral of s_n over the step, which gives
      // the drift and noise_amp weights. The zero mode is white with variance
      // 2T/dt. Cheap, O(n_modes); called whenever the temperature changes.
      //---------------------------------------------------------------------------
      void ou_coefficients(bath_t& b, const double T_scaled){

         const int n_modes = b.n_modes;
         const double dt = b.dt;

         b.ou_decay.resize(n_modes);
         b.ou_diffuse.resize(n_modes);
         b.ou_drift.resize(n_modes);
         b.ou_noise_amp.resize(n_modes);

         for(int n = 0; n < n_modes; n++){
            const double nu = 2.0*M_PI*(n + 1)*T_scaled;
            const double decay = std::exp(-nu*dt);
            const double diffuse = std::sqrt(1.0 - decay*decay);
            b.ou_decay[n] = decay;
            b.ou_diffuse[n] = diffuse;
            b.ou_drift[n] = std::sqrt(2.0*T_scaled/nu)*(decay - 1.0)/dt;
            b.ou_noise_amp[n] = std::sqrt(2.0*T_scaled/nu)*diffuse/dt;
         }

         b.white_amp = std::sqrt(2.0*T_scaled/dt);

         return;

      }

      //---------------------------------------------------------------------------
      // Function to draw one site's sample from the Matsubara bath
      //---------------------------------------------------------------------------
      void ou_draw(bath_t& b, const int site, const double amplitude){

         const int n_modes = b.n_modes;
         const size_t base = static_cast<size_t>(site)*n_modes;

         // unit Gaussians for the modes and the zero mode
         for(int i = 0; i < 3*n_modes + 3; i++) b.scratch[i] = mtrandom::gaussian();

         double out[3];
         ou_update_atom(&b.s_x[base], &b.s_y[base], &b.s_z[base], n_modes,
                        &b.ou_decay[0], &b.ou_diffuse[0], &b.ou_drift[0], &b.ou_noise_amp[0],
                        b.white_amp, &b.scratch[0], out);

         b.x[site] = amplitude*out[0];
         b.y[site] = amplitude*out[1];
         b.z[site] = amplitude*out[2];

         return;

      }

   } // end of internal namespace

} // end of quantum namespace
