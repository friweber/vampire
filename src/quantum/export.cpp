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
#include <fstream>
#include <iostream>

// Vampire headers
#include "material.hpp"
#include "quantum.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      // the export file, opened once on rank 0
      static std::ofstream export_file;

      //---------------------------------------------------------------------------
      // Function to open the export file and write its header. source names
      // the quantity that follows: the auxiliary oscillator of the thermostat
      // or the bath sample added to the effective field.
      //---------------------------------------------------------------------------
      void export_open(const char* source){

         if(vmpi::my_rank != 0) return;

         if(export_atom < 0 || export_atom >= spin_bath.n_sites){
            zlog << zTs() << "Warning: quantum:export-noise-atom " << export_atom
                 << " is not a local atom; noise export disabled." << std::endl;
            export_noise = false;
            return;
         }

         export_file.open(export_filename.c_str(), std::ios::out | std::ios::trunc);
         if(!export_file.is_open()){
            zlog << zTs() << "Warning: cannot open quantum noise export file " << export_filename
                 << "; noise export disabled." << std::endl;
            export_noise = false;
            return;
         }

         export_file << "# VAMPIRE quantum noise export\n"
                     << "# source:    " << source << "\n"
                     << "# spectrum:  " << spectrum_name(spin_bath.spectrum) << "\n"
                     << "# generator: " << generator_name(spin_bath.generator) << "\n"
                     << "# atom:      " << export_atom << "\n"
                     << "# dt:        " << mp::dt_SI << " s\n"
                     << "# time[s]  x  y  z\n";

         return;

      }

      //---------------------------------------------------------------------------
      // Function to append one row: the time at the start of the step and the
      // three components of the exported quantity
      //---------------------------------------------------------------------------
      void export_sample(const double x, const double y, const double z){

         if(vmpi::my_rank != 0 || !export_file.is_open()) return;

         export_file << static_cast<double>(sim::time)*mp::dt_SI << "\t" << x << "\t" << y << "\t" << z << "\n";

         return;

      }

   } // end of internal namespace

} // end of quantum namespace
