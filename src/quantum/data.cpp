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

// Vampire headers
#include "quantum.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   //------------------------------------------------------------------------
   // Shared variables inside quantum module
   //------------------------------------------------------------------------
   namespace internal{

      mode_t mode = inactive;                 // decided in initialize()
      spectrum_t noise_type = classical;      // quantum:noise-type
      bool noise_type_set = false;            // true once the keyword was given
      generator_t noise_generator = on_the_fly; // quantum:noise-generator

      std::vector<internal::mp_t> mp;         // Lorentzian parameters per material

      int n_bath_modes = 30;                  // auxiliary modes per site
      int bath_scan_resolution = 50;          // grid points per axis of the log-bath fit
      int bath_scan_omega_points = 500;       // frequency points of the log-bath fit
      double bath_scan_decades = 3.5;         // half-width of the log-bath search box

      uint64_t window_size = 0;               // 0 = one window spanning the run
      int interpolation_factor = 1;           // fine steps per coarse sample

      bool export_noise = false;              // export the injected noise of one atom
      std::string export_filename = "quantum-noise.dat"; // export file name
      int export_atom = 0;                    // atom to export

      bath_t spin_bath;                       // spin bath
      bath_t lattice_bath;                    // lattice bath (spin-lattice only)

      bool lattice_parameters_set = false;    // spin-lattice handed over its parameters
      std::vector<double> lattice_mass;       // eV ps^2/A^2
      std::vector<double> lattice_damping;    // 1/ps
      std::vector<double> lattice_damping_eq; // 1/ps, during equilibration

      std::vector<double> material_A;         // alpha*omega0^4/gamma
      std::vector<double> material_gamma;     // rad/s
      std::vector<double> material_omega0;    // rad/s
      std::vector<double> material_S0;        // mu_s/(hbar*gamma_e)

      std::vector<double> q_x;                // auxiliary oscillator position
      std::vector<double> q_y;
      std::vector<double> q_z;
      std::vector<double> p_x;                // auxiliary oscillator momentum
      std::vector<double> p_y;
      std::vector<double> p_z;
      std::vector<double> k1;                 // RK4 slopes, [atom*9 + component]
      std::vector<double> k2;
      std::vector<double> k3;
      std::vector<double> k4;
      std::vector<double> y_pred;             // predicted state
      std::vector<double> y_in;               // state at the start of the step

   } // end of internal namespace

} // end of quantum namespace
