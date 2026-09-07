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
//   Input file and material file parameter parsing for the quantum module.
//
//   Supported input file parameters:
//     quantum:noise-type                   classical | quantum | quantum-no-zero
//     quantum:llg-method                   llg-fft | llg-ho
//     quantum:noise-window-size            integer (must be divisible by 6)
//     quantum:noise-interpolation-factor   integer >= 1
//     quantum:bath-modes                   integer >= 2  (HO method, quantum / quantum-no-zero)
//     quantum:bath-scan-resolution         integer >= 5  (quantum-no-zero λ-range scan)
//     quantum:bath-scan-omega-points       integer >= 50 (quantum-no-zero objective ω grid)
//     quantum:bath-scan-decades            double  >= 1.0 (quantum-no-zero search box width)
//     quantum:export-noise                 [filename]
//     quantum:heun-export-noise            [filename]  (direct route only, see below)
//     quantum:heun-export-noise-atom       integer >= 0  (direct route only)
//
//   Supported material file parameters:
//     quantum-lorentzian-width            Gamma [rad/s]
//     quantum-lorentzian-central-frequency  omega0 [rad/s]
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <string>

// Vampire headers
#include "quantum.hpp"
#include "errors.hpp"
#include "vio.hpp"

// Module headers
#include "internal.hpp"
#include "../spinlattice/internal.hpp"   // sld::internal::export_noise* (see heun-export-noise below)

namespace quantum{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for quantum module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="quantum";
      if(key!=prefix) return false;

      // NOTE: enabling the standalone quantum thermostat (internal::enabled) is done
      // per-keyword below, NOT blanket. The bath-tuning keywords (bath-modes,
      // bath-scan-*) are SHARED with the sld module's quantum-noise path and must
      // NOT switch the thermostat on by themselves — otherwise an sld-only run that
      // sets quantum:bath-modes would trigger quantum::initialize() and abort on the
      // missing quantum-lorentzian-* material parameters.

      //------------------------------------------------------------------------
      // Noise spectral density type
      //------------------------------------------------------------------------
      std::string test = "noise-type";
      if(word == test){
         internal::enabled = true;
         test = "classical";
         if(value == test){
            internal::noise_type = internal::classical;
            return true;
         }
         test = "quantum";
         if(value == test){
            internal::noise_type = internal::quantum_zero;
            return true;
         }
         test = "quantum-no-zero";
         if(value == test){
            internal::noise_type = internal::quantum_no_zero;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"classical\"" << std::endl;
            std::cerr << "\t\"quantum\"" << std::endl;
            std::cerr << "\t\"quantum-no-zero\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      //------------------------------------------------------------------------
      // Noise kind for the llg-heun-quantum integrator.
      // Does NOT enable the full quantum thermostat (no Lorentzian params needed).
      //------------------------------------------------------------------------
      test = "heun-noise-type";
      if(word == test){
         test = "classical";
         if(value == test){ internal::heun_noise_kind = quantum::sld_noise::kind_t::classical;          return true; }
         test = "quantum";
         if(value == test){ internal::heun_noise_kind = quantum::sld_noise::kind_t::quantum;            return true; }
         test = "quantum-no-zero";
         if(value == test){ internal::heun_noise_kind = quantum::sld_noise::kind_t::quantum_no_zero;    return true; }
         test = "quantum-fft";
         if(value == test){ internal::heun_noise_kind = quantum::sld_noise::kind_t::quantum_fft;        return true; }
         test = "quantum-no-zero-fft";
         if(value == test){ internal::heun_noise_kind = quantum::sld_noise::kind_t::quantum_no_zero_fft; return true; }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"classical\"" << std::endl;
            std::cerr << "\t\"quantum\"" << std::endl;
            std::cerr << "\t\"quantum-no-zero\"" << std::endl;
            std::cerr << "\t\"quantum-fft\"" << std::endl;
            std::cerr << "\t\"quantum-no-zero-fft\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      //------------------------------------------------------------------------
      // Export noise to file for analysis
      //------------------------------------------------------------------------
      test = "export-noise";
      if(word == test){
         internal::enabled = true;
         internal::export_noise = true;
         // Three accepted forms:
         //   quantum:export-noise                  -> default name
         //   quantum:export-noise = true|1|yes|on  -> default name
         //   quantum:export-noise = mynoise.dat    -> custom name
         // Default name is built in initialize() as <noise-type>_<method>_noise.dat
         // when the filename here is left empty.
         const bool is_bool_true =
            value == "true"  || value == "1"  || value == "yes" ||
            value == "on"    || value == "enable" || value == "enabled";
         if(!value.empty() && !is_bool_true){
            internal::export_noise_filename = value;
         }
         return true;
      }

      //------------------------------------------------------------------------
      // Export noise for the DIRECT route (sim:integrator=llg-heun-quantum).
      //
      // That route reads sld::internal::export_noise (LLGHeun_quantum.cpp),
      // which is normally set by spin-lattice:export-noise -- but ANY
      // "spin-lattice:" keyword unconditionally sets sld::enabled = true
      // (sld/interface.cpp:38), pulling in the full lattice-dynamics module
      // (masses, potentials, phonon setup). For a bare llg-heun-quantum run
      // with no lattice configured that segfaults during SLD's own startup.
      // This keyword pokes the same three flags directly, without touching
      // sld::enabled, so the direct route's noise can be exported on its own.
      //------------------------------------------------------------------------
      test = "heun-export-noise";
      if(word == test){
         // Deliberately does NOT set internal::enabled: that flag also
         // gates quantum::initialize() (initialize_modules.cpp), which the
         // direct route never calls -- turning it on here made the run
         // enter that open-system init path with no omega0/gamma set,
         // which crashed instead of cleanly erroring.
         sld::internal::export_noise = true;
         const bool is_bool_true =
            value == "true"  || value == "1"  || value == "yes" ||
            value == "on"    || value == "enable" || value == "enabled";
         if(!value.empty() && !is_bool_true){
            sld::internal::export_noise_filename = value;
         }
         return true;
      }

      test = "heun-export-noise-atom";
      if(word == test){
         sld::internal::export_noise_atom = vin::str_to_uint64(value);
         return true;
      }

      //------------------------------------------------------------------------
      // LLG integration method selection
      //------------------------------------------------------------------------
      test = "llg-method";
      if(word == test){
         internal::enabled = true;
         test = "llg-ho";
         if(value == test){
            internal::llg_method = internal::llg_ho;
            return true;
         }
         test = "llg-fft";
         if(value == test){
            internal::llg_method = internal::llg_fft;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"llg-ho\"" << std::endl;
            std::cerr << "\t\"llg-fft\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }


      // Specify which parameters is for with integrator method!! (FFT/HO)

      //------------------------------------------------------------------------
      // Noise window size (for windowed FFT noise generation)
      // Must be divisible by OVERLAP_SAVE_SEGMENTS for the overlap-save scheme.
      //------------------------------------------------------------------------
      test = "noise-window-size";
      if(word == test){
         internal::enabled = true;
         uint64_t ws = vin::str_to_uint64(value);
         vin::check_for_valid_int(ws, word, line, prefix, uint64_t(1), uint64_t(6000000000), "input", "> 0");
         if(ws % internal::OVERLAP_SAVE_SEGMENTS != 0){
            std::cerr << "Error: Quantum window size must be divisible by "
                      << internal::OVERLAP_SAVE_SEGMENTS
                      << " (overlap-save segment count)." << std::endl;
            return false;
         }
         internal::window_size = ws;
         return true;
      }

      //------------------------------------------------------------------------
      // Interpolation factor M (coarse-to-fine time step ratio)
      //------------------------------------------------------------------------
      test = "noise-interpolation-factor";
      if(word == test){
         internal::enabled = true;
         int md = vin::str_to_uint64(value);
         vin::check_for_valid_int(md, word, line, prefix, 1, 1000000, "input", "> 0");
         internal::M_decimation = md;
         return true;
      }

      //------------------------------------------------------------------------
      // Number of auxiliary bath modes per atom (HO method).
      // Used by both quantum (orn-uhl, Matsubara modes) and
      // quantum-no-zero (log-bath-opt-range) paths.
      //------------------------------------------------------------------------
      test = "bath-modes";
      if(word == test){
         int n = vin::str_to_uint64(value);
         vin::check_for_valid_int(n, word, line, prefix, 2, 100000, "input", "> 1");
         internal::n_bath_modes = n;
         return true;
      }

      //------------------------------------------------------------------------
      // Grid points per axis for the (log_ls, log_le) λ-range scan
      // used in setup_log_bath_opt (quantum-no-zero noise type).
      // Default 50 reproduces the cmp_noise reference fit quality.
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // Depth of the cascade OU chain per log-bath mode (quantum-no-zero only).
      // M=1 = no cascade (default, current behaviour).
      // M>1: each mode drives a chain of M coupled OU processes; output from
      // the last level gives ω^{-2M} high-frequency roll-off.
      //------------------------------------------------------------------------
      test = "cascade-modes";
      if(word == test){
         int m = vin::str_to_uint64(value);
         vin::check_for_valid_int(m, word, line, prefix, 1, 100, "input", "in [1, 100]");
         internal::n_cascade_modes = m;
         return true;
      }

      //------------------------------------------------------------------------
      // LP cutoff for the Butterworth post-filter (quantum-no-zero only).
      // When set, quantum_no_zero uses generate_quantum_noise_log_bath_opt_filtered
      // instead of the raw log-bath or cascade generator.
      // Accepts a floating-point value with optional !THz unit suffix.
      // A value of 0 (default) disables the filter entirely.
      //------------------------------------------------------------------------
      test = "filter-cutoff";
      if(word == test){
         double fc = vin::str_to_double(value);
         vin::check_for_valid_value(fc, word, line, prefix, unit, "frequency",
                                    0.0, 1.0e15, "input", ">= 0");
         internal::butter_cutoff_Hz    = fc;   // stored in Hz after unit conversion
         internal::butter_coeffs_valid = false; // force re-compute on next generate call
         return true;
      }

      test = "bath-scan-resolution";
      if(word == test){
         int n = vin::str_to_uint64(value);
         vin::check_for_valid_int(n, word, line, prefix, 5, 500, "input", "in [5, 500]");
         internal::bath_scan_resolution = n;
         return true;
      }

      //------------------------------------------------------------------------
      // Number of points in the linear ω grid used to evaluate the SSE
      // objective in setup_log_bath_opt (quantum-no-zero noise type).
      // Default 500 matches cmp_noise; raise for finer fits, lower for
      // faster initialisation.
      //------------------------------------------------------------------------
      test = "bath-scan-omega-points";
      if(word == test){
         int n = vin::str_to_uint64(value);
         vin::check_for_valid_int(n, word, line, prefix, 50, 5000, "input", "in [50, 5000]");
         internal::bath_scan_omega_points = n;
         return true;
      }

      //------------------------------------------------------------------------
      // Width (in log-decades) of the asymmetric search box around log10(T)
      // and log10(max(T, ω₀)) used in setup_log_bath_opt (quantum-no-zero).
      // Larger -> explore further into small-λ and large-λ regimes.
      //------------------------------------------------------------------------
      test = "bath-scan-decades";
      if(word == test){
         double d = vin::str_to_double(value);
         vin::check_for_valid_value(d, word, line, prefix, unit, "none", 1.0, 10.0, "input", "in [1.0, 10.0]");
         internal::bath_scan_decades = d;
         return true;
      }

      return false;
   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      std::string prefix="material:";

      // Resize mp vector if necessary
      if(internal::mp.size() <= (unsigned int)super_index){
         internal::mp.resize(super_index + 1);
      }

      //------------------------------------------------------------------------
      // Lorentzian width Gamma [rad/s]
      //------------------------------------------------------------------------
      std::string test = "quantum-lorentzian-width";
      if(word == test){
         internal::enabled = true;
         double gamma = vin::str_to_double(value);
         vin::check_for_valid_value(gamma, word, line, prefix, unit, "frequency", 0.0, 1.0e15, "material", "> 0");
         internal::mp[super_index].gamma.set(gamma);
         return true;
      }

      //------------------------------------------------------------------------
      // Lorentzian central frequency omega0 [rad/s]
      //------------------------------------------------------------------------
      test = "quantum-lorentzian-central-frequency";
      if(word == test){
         internal::enabled = true;
         double omega0 = vin::str_to_double(value);
         vin::check_for_valid_value(omega0, word, line, prefix, unit, "frequency", 0.0, 1.0e15, "material", "> 0");
         internal::mp[super_index].omega0.set(omega0);
         return true;
      }

      return false;
   }

} // end of quantum namespace
