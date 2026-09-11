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

#ifndef QUANTUM_INTERNAL_H_
#define QUANTUM_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the quantum module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>
#include <string>
#include <vector>

#ifdef FFT
#include <fftw3.h>
#endif

// Vampire headers
#include "quantum.hpp"

namespace quantum{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      // What the module does in this run, decided once in initialize()
      enum mode_t{
         inactive,   // integrator uses its own white thermal noise
         thermostat, // open-system LLG with auxiliary oscillator (llg-quantum)
         bath        // coloured noise handed to llg-heun or spin-lattice
      };

      // Spectral density of the bath, S(w) in units of the reduced temperature
      enum spectrum_t{
         classical,       // 2T
         quantum_zero,    // w coth(w/2T)        (with zero-point term)
         quantum_no_zero  // w (coth(w/2T) - 1)  (without zero-point term)
      };

      // How the coloured samples are produced
      enum generator_t{
         on_the_fly,    // auxiliary Ornstein-Uhlenbeck modes, one update per step
         pre_generated  // spectrally shaped white noise, FFT in windows
      };

      // simple initialised class for set variables
      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifying variable has been set

      public:
         set_double_t() : value(0.0), setf(false) { }
         void set(double in_value){ value = in_value; setf = true; };
         double get(){ return value; };
         bool is_set(){ return setf; };
      };

      // per-material parameters of the thermostat's Lorentzian response
      class mp_t{

         private:

         public:
            set_double_t gamma;  // Lorentzian width (rad/s)
            set_double_t omega0; // Lorentzian central frequency (rad/s)

            // constructor
            mp_t (){
               gamma.set(0.0);
               omega0.set(0.0);
            };
      };

      // Window state of the pre-generated (overlap-save FFT) generator. The
      // FFT length is six segments; the first and last segment of every
      // transform are discarded, so a window yields four segments of
      // usable coarse samples. coarse[] holds one realisation per site and
      // component, laid out as [(3*site + component)*n_valid + j].
      class fft_t{

         public:
            int M;                 // fine steps per coarse sample
            int n_coarse;          // FFT length, 6*seg
            int seg;               // segment length in coarse samples
            int n_valid;           // 4*seg coarse samples kept per window
            bool primed;           // false until the first window exists
            uint64_t step;         // fine step counter of this bath
            uint64_t window_start; // fine step at which the current window starts
            double norm;           // output scale, 1/(n_coarse*sqrt(M*dt))
            std::vector< std::vector<double> > sqrt_psd; // [material][n_coarse/2+1]
            std::vector<double> tail;   // 2*seg white samples per realisation, carried over
            std::vector<double> coarse; // shaped coarse samples of the current window
            #ifdef FFT
            fftw_plan forward;
            fftw_plan backward;
            double* in;
            fftw_complex* out;
            double* result;
            #endif

            // constructor
            fft_t () : M(1), n_coarse(0), seg(0), n_valid(0), primed(false),
                       step(0), window_start(0), norm(0.0)
            {
               #ifdef FFT
               forward = 0;
               backward = 0;
               in = 0;
               out = 0;
               result = 0;
               #endif
            };
      };

      // One coloured bath: n_sites sites, three components each. The same
      // type serves the spin bath of every consumer and the lattice bath of
      // the spin-lattice module; only dt, the temperature scale and the
      // amplitudes differ.
      class bath_t{

         public:
            spectrum_t spectrum;
            generator_t generator;
            bool lorentzian;       // shape the pre-generated spectrum by the thermostat Lorentzian
            int n_sites;
            int n_modes;           // auxiliary modes per site (on-the-fly generator)
            double dt;             // time step in the bath's own units
            double T_scale;        // Kelvin to reduced temperature
            double last_T_scaled;  // temperature the coefficients were built for
            std::vector<double> amp;    // per-material output amplitude
            std::vector<double> amp_eq; // per-material output amplitude during equilibration
            // Matsubara Ornstein-Uhlenbeck modes (quantum_zero)
            std::vector<double> ou_decay;
            std::vector<double> ou_diffuse;
            std::vector<double> ou_drift;
            std::vector<double> ou_noise_amp;
            double white_amp;
            // log-spaced Ornstein-Uhlenbeck bath (quantum_no_zero)
            std::vector<double> lb_lambda;
            std::vector<double> lb_decay;
            std::vector<double> lb_coeff;
            std::vector<double> lb_amp;
            // per-site mode state, [site*n_modes + mode], shared by both kinds
            std::vector<double> s_x;
            std::vector<double> s_y;
            std::vector<double> s_z;
            std::vector<double> scratch; // Gaussian draws for one site
            fft_t fft;
            // this step's samples
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> z;

            // constructor
            bath_t () : spectrum(classical), generator(on_the_fly), lorentzian(false),
                        n_sites(0), n_modes(0), dt(0.0), T_scale(0.0),
                        last_T_scaled(0.0), white_amp(0.0) { };
      };

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern mode_t mode;                 // what the module does in this run
      extern spectrum_t noise_type;       // requested spectrum (quantum:noise-type)
      extern bool noise_type_set;         // true once quantum:noise-type was given
      extern generator_t noise_generator; // requested generator (quantum:noise-generator)

      extern std::vector<internal::mp_t> mp; // Lorentzian parameters (thermostat only)

      extern int n_bath_modes;            // auxiliary modes per site
      extern int bath_scan_resolution;    // grid points per axis of the log-bath fit
      extern int bath_scan_omega_points;  // frequency points of the log-bath fit
      extern double bath_scan_decades;    // width of the log-bath search box (decades)

      extern uint64_t window_size;        // pre-generated window in fine steps (0 = whole run)
      extern int interpolation_factor;    // fine steps per coarse sample

      extern bool export_noise;           // write the injected noise of one atom
      extern std::string export_filename; // file for the noise export
      extern int export_atom;             // atom whose noise is exported

      extern bath_t spin_bath;            // spin bath (thermostat, llg-heun, spin-lattice)
      extern bath_t lattice_bath;         // lattice bath (spin-lattice only)

      extern bool lattice_parameters_set;            // set_lattice_parameters() was called
      extern std::vector<double> lattice_mass;       // per-material atomic mass (eV ps^2/A^2)
      extern std::vector<double> lattice_damping;    // per-material lattice damping (1/ps)
      extern std::vector<double> lattice_damping_eq; // same during equilibration

      // thermostat: per-material Lorentzian data
      extern std::vector<double> material_A;      // alpha*omega0^4/gamma
      extern std::vector<double> material_gamma;  // Lorentzian width (rad/s)
      extern std::vector<double> material_omega0; // Lorentzian central frequency (rad/s)
      extern std::vector<double> material_S0;     // moment in units of hbar*gamma_e

      // thermostat: per-atom auxiliary oscillator and RK4 storage
      extern std::vector<double> q_x;
      extern std::vector<double> q_y;
      extern std::vector<double> q_z;
      extern std::vector<double> p_x;
      extern std::vector<double> p_y;
      extern std::vector<double> p_z;
      extern std::vector<double> k1;     // [atom*9 + component]
      extern std::vector<double> k2;
      extern std::vector<double> k3;
      extern std::vector<double> k4;
      extern std::vector<double> y_pred;
      extern std::vector<double> y_in;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      // initialize.cpp
      int num_local_atoms();
      bool run_length(uint64_t& n_steps);
      bool dynamic_temperature_program();
      const char* spectrum_name(const spectrum_t s);
      const char* generator_name(const generator_t g);

      // bath.cpp
      void setup_bath(bath_t& b, const spectrum_t spectrum, const generator_t generator,
                      const int n_sites, const double dt, const double T_scale,
                      const std::vector<double>& amp, const std::vector<double>& amp_eq,
                      const bool lorentzian, const uint64_t n_run);
      void refresh_temperature(bath_t& b);
      void draw(bath_t& b);
      void draw_white_thermostat(bath_t& b);

      // bath_ou.cpp
      void ou_coefficients(bath_t& b, const double T_scaled);
      void ou_draw(bath_t& b, const int site, const double amplitude);

      // bath_log.cpp
      void log_bath_fit(bath_t& b, const double T_scaled);
      void log_bath_coefficients(bath_t& b, const double T_scaled);
      void log_bath_draw(bath_t& b, const int site, const double amplitude);

      // bath_fft.cpp
      void fft_setup(bath_t& b, const uint64_t n_run);
      void fft_generate_window(bath_t& b);
      void fft_draw(bath_t& b, const std::vector<double>& amplitude);
      void fft_release(bath_t& b);

      // export.cpp
      void export_open(const char* source);
      void export_sample(const double x, const double y, const double z);

      // llg.cpp
      void allocate_thermostat(const int n_atoms);
      void rk4_save(const int start_index, const int end_index);
      void rk4_stage(const int start_index, const int end_index, const int stage, const bool ho);
      void rk4_writeback_spin(const int start_index, const int end_index);
      void rk4_finish(const int start_index, const int end_index);
      void export_thermostat_sample();

      // llg_serial.cpp, llg_mpi.cpp
      void llg_ho_serial();
      void llg_fft_serial();
      void llg_ho_mpi();
      void llg_fft_mpi();

   } // end of internal namespace

} // end of quantum namespace

#endif //QUANTUM_INTERNAL_H_
