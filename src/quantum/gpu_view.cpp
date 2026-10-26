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

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "gpu.hpp"
#include "material.hpp"
#include "sim.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   //---------------------------------------------------------------------------
   // Function to report whether the GPU takes over the module's work
   //---------------------------------------------------------------------------
   bool gpu_active(){
      return gpu::acceleration && internal::mode != internal::inactive;
   }

   bool gpu_bath(){
      return gpu::acceleration && internal::mode == internal::bath;
   }

   bool gpu_thermostat(){
      return gpu::acceleration && internal::mode == internal::thermostat;
   }

   //---------------------------------------------------------------------------
   // Function to fill the GPU view for this step. The coefficients are rebuilt
   // when sim::temperature has moved (exactly what the CPU draw does), and the
   // per-material amplitude follows the equilibration / run phase.
   //---------------------------------------------------------------------------
   void gpu_view(gpu_view_t& v){

      using namespace internal;

      bath_t& b = spin_bath;

      v.mode = (mode == thermostat) ? 1 : (mode == bath) ? 2 : 0;
      v.spectrum = (b.spectrum == classical) ? 0 : (b.spectrum == quantum_zero) ? 1 : 2;
      v.n_modes = b.n_modes;
      v.n_materials = static_cast<int>(mp::material.size());
      v.dt = mp::dt;

      if(b.spectrum != classical) refresh_temperature(b);
      v.T_scaled = b.last_T_scaled;

      if(b.spectrum == quantum_zero){
         v.c0 = b.ou_decay.data();
         v.c1 = b.ou_diffuse.data();
         v.c2 = b.ou_drift.data();
         v.c3 = b.ou_noise_amp.data();
         v.white_amp = b.white_amp;
      }
      else{
         v.c0 = b.lb_decay.data();
         v.c1 = b.lb_coeff.data();
         v.c2 = b.lb_amp.data();
         v.c3 = 0;
         v.white_amp = 0.0;
      }

      // equilibration uses its own amplitudes, as in draw()
      v.amplitude = (sim::time < sim::equilibration_time) ? b.amp_eq.data() : b.amp.data();

      v.A = material_A.data();
      v.gamma = material_gamma.data();
      v.omega0 = material_omega0.data();

      // combined thermostat: the Ohmic bath has the spectrum, temperature and step of the
      // spin bath, hence the same coefficients; only its amplitudes differ
      v.ohmic = (mode == thermostat && ohmic) ? 1 : 0;
      v.alpha_e = material_alpha_e.data();
      v.amplitude_e = ohmic_bath.amp.data();

      // classical open-system thermostat: white noise with the prefactor of draw_white_thermostat()
      v.white.assign(v.n_materials, 0.0);
      v.white_e.assign(v.n_materials, 0.0);
      if(mode == thermostat && b.spectrum == classical){
         for(int m = 0; m < v.n_materials; m++){
            double T = sim::temperature;
            if(sim::local_temperature) T = mp::material[m].temperature;
            const double a  = mp::material[m].temperature_rescaling_alpha;
            const double Tc = mp::material[m].temperature_rescaling_Tc;
            const double T_rescaled = (T < Tc) ? Tc*std::pow(T/Tc, a) : T;
            const double T_scaled = T_rescaled*b.T_scale;
            v.white[m] = std::sqrt(2.0*material_gamma[m]*material_A[m]*T_scaled/(material_S0[m]*b.dt));
            if(v.ohmic) v.white_e[m] = std::sqrt(2.0*material_alpha_e[m]*T_scaled/(material_S0[m]*b.dt));
         }
      }

      return;

   }

} // end of quantum namespace
