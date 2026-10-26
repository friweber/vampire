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

#ifndef QUANTUM_LLG_ATOM_H_
#define QUANTUM_LLG_ATOM_H_
//
//---------------------------------------------------------------------
// Per-atom arithmetic of the quantum module: the equations of motion,
// the RK4 stage updates and the auxiliary-bath updates for one site.
// Everything here works on raw pointers and scalars only, so the same
// functions serve the CPU steppers and a GPU kernel. Do not add
// std::vector, mtrandom or sim:: here.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cmath>

#ifdef __CUDACC__
#define QUANTUM_HOST_DEVICE __host__ __device__
#else
#define QUANTUM_HOST_DEVICE
#endif

namespace quantum{

   namespace internal{

      //-------------------------------------------------------------------------
      // Equation of motion with the bath sample driving the oscillator.
      // State y = (S, q, p); the spin precesses in H + q, the oscillator is
      // driven by the spin and by the (unfiltered) noise:
      //    dS/dt = S x (H + q)
      //    dq/dt = p
      //    dp/dt = -omega0^2 q - Gamma p + A S + noise
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void ho_equation_of_motion(const double* y, const double* H,
                                                            const double* noise,
                                                            const double A, const double Gamma,
                                                            const double omega0, double* dydt){

         dydt[0] = y[1]*(H[2] + y[5]) - y[2]*(H[1] + y[4]);
         dydt[1] = y[2]*(H[0] + y[3]) - y[0]*(H[2] + y[5]);
         dydt[2] = y[0]*(H[1] + y[4]) - y[1]*(H[0] + y[3]);

         dydt[3] = y[6];
         dydt[4] = y[7];
         dydt[5] = y[8];

         dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0] + noise[0];
         dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1] + noise[1];
         dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2] + noise[2];

         return;

      }

      //-------------------------------------------------------------------------
      // Equation of motion of the combined thermostat: the Lorentzian
      // oscillator above (phonon bath) plus an Ohmic bath acting on the spin
      // directly (electron bath), a Gilbert damping alpha_e with its own
      // noise field b. In Gilbert form, with H' = H + q + b,
      //    dS/dt = S x H' - alpha_e S x dS/dt
      // solved for dS/dt (Landau-Lifshitz form):
      //    dS/dt = [S x H' - alpha_e S x (S x H')]/(1 + alpha_e^2)
      // The oscillator is unchanged. At low frequency q contains
      // -alpha dS/dt, so the total Gilbert damping there is alpha + alpha_e.
      // alpha_e = 0, b = 0 is ho_equation_of_motion.
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void ho_ohmic_equation_of_motion(const double* y, const double* H,
                                                                  const double* noise, const double* b,
                                                                  const double A, const double Gamma,
                                                                  const double omega0, const double alpha_e,
                                                                  double* dydt){

         const double hx = H[0] + y[3] + b[0];
         const double hy = H[1] + y[4] + b[1];
         const double hz = H[2] + y[5] + b[2];

         // S x H'
         const double tx = y[1]*hz - y[2]*hy;
         const double ty = y[2]*hx - y[0]*hz;
         const double tz = y[0]*hy - y[1]*hx;

         // S x (S x H')
         const double dx = y[1]*tz - y[2]*ty;
         const double dy = y[2]*tx - y[0]*tz;
         const double dz = y[0]*ty - y[1]*tx;

         const double f = 1.0/(1.0 + alpha_e*alpha_e);
         dydt[0] = f*(tx - alpha_e*dx);
         dydt[1] = f*(ty - alpha_e*dy);
         dydt[2] = f*(tz - alpha_e*dz);

         dydt[3] = y[6];
         dydt[4] = y[7];
         dydt[5] = y[8];

         dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0] + noise[0];
         dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1] + noise[1];
         dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2] + noise[2];

         return;

      }

      //-------------------------------------------------------------------------
      // Equation of motion for pre-generated noise. The sample is already
      // shaped by the Lorentzian and enters through H; the oscillator runs
      // without a noise term.
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void fft_equation_of_motion(const double* y, const double* H,
                                                             const double A, const double Gamma,
                                                             const double omega0, double* dydt){

         dydt[0] = y[1]*(H[2] + y[5]) - y[2]*(H[1] + y[4]);
         dydt[1] = y[2]*(H[0] + y[3]) - y[0]*(H[2] + y[5]);
         dydt[2] = y[0]*(H[1] + y[4]) - y[1]*(H[0] + y[3]);

         dydt[3] = y[6];
         dydt[4] = y[7];
         dydt[5] = y[8];

         dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0];
         dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1];
         dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2];

         return;

      }

      //-------------------------------------------------------------------------
      // Normalise the spin part of a 9-component state in place
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void renormalise_spin(double* y){

         const double inv_mag = 1.0/sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
         y[0] *= inv_mag;
         y[1] *= inv_mag;
         y[2] *= inv_mag;

         return;

      }

      //-------------------------------------------------------------------------
      // RK4 predictor: y_pred = y_in + coeff*k, spin renormalised
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void rk4_predict(const double* y_in, const double* k,
                                                  const double coeff, double* y_pred){

         for(int i = 0; i < 9; i++) y_pred[i] = y_in[i] + coeff*k[i];
         renormalise_spin(y_pred);

         return;

      }

      //-------------------------------------------------------------------------
      // RK4 combination: y_out = y_in + dt/6 (k1 + 2 k2 + 2 k3 + k4), spin renormalised
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void rk4_combine(const double* y_in,
                                                  const double* k1, const double* k2,
                                                  const double* k3, const double* k4,
                                                  const double dt_over_6, double* y_out){

         for(int i = 0; i < 9; i++){
            y_out[i] = y_in[i] + dt_over_6*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
         }
         renormalise_spin(y_out);

         return;

      }

      //-------------------------------------------------------------------------
      // One step of the Matsubara Ornstein-Uhlenbeck bath for one site.
      // The modes s_n are updated exactly, s_n <- decay_n s_n + diffuse_n xi_n,
      // and the noise delivered over the step is the integral of the modes
      // plus a white zero-mode term:
      //    out = white_amp xi_0 + sum_n (drift_n s_n(old) + noise_amp_n xi_n)
      // xi holds 3*n_modes + 3 unit Gaussians: x modes, y modes, z modes,
      // then the three zero-mode draws.
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void ou_update_atom(double* sx, double* sy, double* sz,
                                                     const int n_modes,
                                                     const double* decay, const double* diffuse,
                                                     const double* drift, const double* noise_amp,
                                                     const double white_amp,
                                                     const double* xi, double* out){

         double aux_x = 0.0;
         double aux_y = 0.0;
         double aux_z = 0.0;

         for(int n = 0; n < n_modes; n++){

            const double zx = xi[n];
            const double zy = xi[n_modes + n];
            const double zz = xi[2*n_modes + n];

            // contribution of this mode, using the state before the update
            aux_x += drift[n]*sx[n] + noise_amp[n]*zx;
            aux_y += drift[n]*sy[n] + noise_amp[n]*zy;
            aux_z += drift[n]*sz[n] + noise_amp[n]*zz;

            // exact update of the mode
            sx[n] = decay[n]*sx[n] + diffuse[n]*zx;
            sy[n] = decay[n]*sy[n] + diffuse[n]*zy;
            sz[n] = decay[n]*sz[n] + diffuse[n]*zz;

         }

         out[0] = white_amp*xi[3*n_modes + 0] + aux_x;
         out[1] = white_amp*xi[3*n_modes + 1] + aux_y;
         out[2] = white_amp*xi[3*n_modes + 2] + aux_z;

         return;

      }

      //-------------------------------------------------------------------------
      // One step of the log-spaced Ornstein-Uhlenbeck bath for one site.
      // Each mode s_k relaxes with rate lambda_k, s_k <- decay_k s_k + coeff_k xi_k,
      // and the noise is the weighted sum of the updated modes:
      //    out = sum_k amp_k s_k(new)
      // xi holds 3*n_modes unit Gaussians: x modes, y modes, z modes.
      //-------------------------------------------------------------------------
      QUANTUM_HOST_DEVICE inline void log_bath_update_atom(double* sx, double* sy, double* sz,
                                                           const int n_modes,
                                                           const double* decay, const double* coeff,
                                                           const double* amp,
                                                           const double* xi, double* out){

         double aux_x = 0.0;
         double aux_y = 0.0;
         double aux_z = 0.0;

         for(int k = 0; k < n_modes; k++){

            sx[k] = decay[k]*sx[k] + coeff[k]*xi[k];
            sy[k] = decay[k]*sy[k] + coeff[k]*xi[n_modes + k];
            sz[k] = decay[k]*sz[k] + coeff[k]*xi[2*n_modes + k];

            aux_x += amp[k]*sx[k];
            aux_y += amp[k]*sy[k];
            aux_z += amp[k]*sz[k];

         }

         out[0] = aux_x;
         out[1] = aux_y;
         out[2] = aux_z;

         return;

      }

   } // end of internal namespace

} // end of quantum namespace

#endif //QUANTUM_LLG_ATOM_H_
