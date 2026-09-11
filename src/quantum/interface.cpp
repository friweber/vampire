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
#include <string>

// Vampire headers
#include "errors.hpp"
#include "quantum.hpp"
#include "vio.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      // note a deprecated spelling once in the log
      static void deprecated(const std::string& old_key, const std::string& new_key){
         zlog << zTs() << "Warning: " << old_key << " is deprecated, use " << new_key << std::endl;
         return;
      }

      // true for the bare forms of a switch keyword
      static bool is_true(const std::string& value){
         return value.empty() || value == "true" || value == "1" || value == "yes" ||
                value == "on" || value == "enable" || value == "enabled";
      }

      //---------------------------------------------------------------------------
      // Function to set the spectrum from its name. The old *-fft values are
      // accepted and select the pre-generated generator as well.
      //---------------------------------------------------------------------------
      static bool set_noise_type(const std::string& value, const std::string& word){

         noise_type_set = true;

         std::string test = "classical";
         if(value == test){ noise_type = classical; return true; }
         test = "quantum";
         if(value == test){ noise_type = quantum_zero; return true; }
         test = "quantum-no-zero";
         if(value == test){ noise_type = quantum_no_zero; return true; }
         test = "quantum-fft";
         if(value == test){
            deprecated("quantum:noise-type = quantum-fft", "quantum:noise-type = quantum with quantum:noise-generator = pre-generated");
            noise_type = quantum_zero;
            noise_generator = pre_generated;
            return true;
         }
         test = "quantum-no-zero-fft";
         if(value == test){
            deprecated("quantum:noise-type = quantum-no-zero-fft", "quantum:noise-type = quantum-no-zero with quantum:noise-generator = pre-generated");
            noise_type = quantum_no_zero;
            noise_generator = pre_generated;
            return true;
         }

         terminaltextcolor(RED);
         std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"classical\"" << std::endl;
         std::cerr << "\t\"quantum\"" << std::endl;
         std::cerr << "\t\"quantum-no-zero\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
         return false;

      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Function to process input file parameters for quantum module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="quantum";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="noise-type";
      if(word==test){
         return internal::set_noise_type(value, word);
      }
      //--------------------------------------------------------------------
      test="noise-generator";
      if(word==test){
         test="on-the-fly";
         if(value==test){
            internal::noise_generator = internal::on_the_fly;
            return true;
         }
         test="pre-generated";
         if(value==test){
            internal::noise_generator = internal::pre_generated;
            return true;
         }
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"on-the-fly\"" << std::endl;
         std::cerr << "\t\"pre-generated\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      //--------------------------------------------------------------------
      test="noise-window-size";
      if(word==test){
         uint64_t ws = vin::str_to_uint64(value);
         vin::check_for_valid_int(ws, word, line, prefix, uint64_t(0), uint64_t(6000000000), "input", "0 (whole run) or > 0");
         internal::window_size = ws;
         return true;
      }
      //--------------------------------------------------------------------
      test="noise-interpolation-factor";
      if(word==test){
         int md = vin::str_to_uint64(value);
         vin::check_for_valid_int(md, word, line, prefix, 1, 1000000, "input", "> 0");
         internal::interpolation_factor = md;
         return true;
      }
      //--------------------------------------------------------------------
      test="bath-modes";
      if(word==test){
         int n = vin::str_to_uint64(value);
         vin::check_for_valid_int(n, word, line, prefix, 2, 100000, "input", "> 1");
         internal::n_bath_modes = n;
         return true;
      }
      //--------------------------------------------------------------------
      test="bath-scan-resolution";
      if(word==test){
         int n = vin::str_to_uint64(value);
         vin::check_for_valid_int(n, word, line, prefix, 5, 500, "input", "in [5, 500]");
         internal::bath_scan_resolution = n;
         return true;
      }
      //--------------------------------------------------------------------
      test="bath-scan-omega-points";
      if(word==test){
         int n = vin::str_to_uint64(value);
         vin::check_for_valid_int(n, word, line, prefix, 50, 5000, "input", "in [50, 5000]");
         internal::bath_scan_omega_points = n;
         return true;
      }
      //--------------------------------------------------------------------
      test="bath-scan-decades";
      if(word==test){
         double d = vin::str_to_double(value);
         vin::check_for_valid_value(d, word, line, prefix, unit, "none", 1.0, 10.0, "input", "in [1.0, 10.0]");
         internal::bath_scan_decades = d;
         return true;
      }
      //--------------------------------------------------------------------
      test="export-noise";
      if(word==test){
         internal::export_noise = true;
         if(!internal::is_true(value)) internal::export_filename = value;
         return true;
      }
      //--------------------------------------------------------------------
      test="export-noise-atom";
      if(word==test){
         int a = vin::str_to_uint64(value);
         vin::check_for_valid_int(a, word, line, prefix, 0, 100000000, "input", ">= 0");
         internal::export_atom = a;
         return true;
      }
      //--------------------------------------------------------------------
      // Deprecated spellings, kept so older input files keep running
      //--------------------------------------------------------------------
      test="heun-noise-type";
      if(word==test){
         internal::deprecated("quantum:heun-noise-type", "quantum:noise-type");
         return internal::set_noise_type(value, word);
      }
      //--------------------------------------------------------------------
      test="llg-method";
      if(word==test){
         internal::deprecated("quantum:llg-method", "quantum:noise-generator");
         test="llg-ho";
         if(value==test){
            internal::noise_generator = internal::on_the_fly;
            return true;
         }
         test="llg-fft";
         if(value==test){
            internal::noise_generator = internal::pre_generated;
            return true;
         }
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"llg-ho\"" << std::endl;
         std::cerr << "\t\"llg-fft\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      //--------------------------------------------------------------------
      test="heun-export-noise";
      if(word==test){
         internal::deprecated("quantum:heun-export-noise", "quantum:export-noise");
         internal::export_noise = true;
         if(!internal::is_true(value)) internal::export_filename = value;
         return true;
      }
      //--------------------------------------------------------------------
      test="heun-export-noise-atom";
      if(word==test){
         internal::deprecated("quantum:heun-export-noise-atom", "quantum:export-noise-atom");
         int a = vin::str_to_uint64(value);
         vin::check_for_valid_int(a, word, line, prefix, 0, 100000000, "input", ">= 0");
         internal::export_atom = a;
         return true;
      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      std::string prefix="material:";

      // Resize mp vector if necessary
      if(internal::mp.size() <= static_cast<unsigned int>(super_index)){
         internal::mp.resize(super_index + 1);
      }

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="quantum-lorentzian-width";
      if(word==test){
         double gamma = vin::str_to_double(value);
         vin::check_for_valid_value(gamma, word, line, prefix, unit, "none", 0.0, 1.0e15, "material", "> 0");
         internal::mp[super_index].gamma.set(gamma);
         return true;
      }
      //--------------------------------------------------------------------
      test="quantum-lorentzian-central-frequency";
      if(word==test){
         double omega0 = vin::str_to_double(value);
         vin::check_for_valid_value(omega0, word, line, prefix, unit, "none", 0.0, 1.0e15, "material", "> 0");
         internal::mp[super_index].omega0.set(omega0);
         return true;
      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of quantum namespace
