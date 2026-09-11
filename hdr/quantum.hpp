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

#ifndef QUANTUM_H_
#define QUANTUM_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "quantum.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for the quantum (coloured noise) module
//--------------------------------------------------------------------------------
namespace quantum{

   //-----------------------------------------------------------------------------
   // Function to initialise the quantum module. Decides from sim::integrator
   // and quantum:noise-type whether the module runs the open-system
   // thermostat, supplies a coloured bath to another integrator, or stays
   // inactive.
   //-----------------------------------------------------------------------------
   void initialize();

   //-----------------------------------------------------------------------------
   // Function to release FFTW plans and buffers at program exit
   //-----------------------------------------------------------------------------
   void cleanup();

   //-----------------------------------------------------------------------------
   // Function to integrate one open-system LLG step (sim:integrator = llg-quantum)
   //-----------------------------------------------------------------------------
   void llg();

   //-----------------------------------------------------------------------------
   // Coloured spin bath for other integrators (llg-heun, spin-lattice).
   // enabled() is true when the bath replaces the classical thermal field;
   // the caller then draws once per step with generate() and adds field()
   // to its effective field.
   //-----------------------------------------------------------------------------
   bool enabled();
   void generate();
   double field(const int atom, const int component);
   void add_field(const int start_index, const int end_index);

   //-----------------------------------------------------------------------------
   // Lattice (phonon) bath for the spin-lattice module. The lattice module
   // hands over its per-material parameters before quantum::initialize()
   // runs; lattice_field() returns this step's force noise for one atom.
   //-----------------------------------------------------------------------------
   void set_lattice_parameters(const std::vector<double>& mass,
                               const std::vector<double>& damping,
                               const std::vector<double>& damping_eq);
   double lattice_field(const int atom, const int component);

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for quantum module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word,
                              std::string const value, std::string const unit,
                              int const line);

   //-----------------------------------------------------------------------------
   // Function to process material parameters
   //-----------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value,
                                 std::string const unit, int const line,
                                 int const super_index, const int sub_index);

} // end of quantum namespace

#endif //QUANTUM_H_
