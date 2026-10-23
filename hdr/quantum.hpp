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
   // enabled() is true when the bath replaces the classical thermal field. The
   // integrator draws once per step with generate(); llg-heun then receives the
   // sample through sim::calculate_thermal_fields(), which calls
   // thermal_fields(), and spin-lattice reads it per atom with field().
   //-----------------------------------------------------------------------------
   bool enabled();
   void generate();
   double field(const int atom, const int component);
   void thermal_fields(const int start_index, const int end_index);

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

   //---------------------------------------------------------------------------
   // GPU interface. The CUDA module integrates the coloured bath and the
   // open-system LLG on the device; the host keeps the bath coefficients
   // (built for the current temperature) and hands them over through this
   // view once per step. Nothing here draws random numbers.
   //---------------------------------------------------------------------------
   struct gpu_view_t{
      int mode;            // 0 inactive, 1 open-system LLG (llg-quantum), 2 coloured bath (llg-heun)
      int spectrum;        // 0 classical, 1 quantum, 2 quantum-no-zero
      int n_modes;         // auxiliary modes per site and component
      int n_materials;
      double T_scaled;     // reduced temperature the coefficients belong to
      double dt;           // reduced time step (mp::dt)
      double white_amp;    // Matsubara bath: amplitude of the white zero mode
      const double* c0;    // Matsubara: decay      | log bath: decay
      const double* c1;    //            diffuse    |           coeff
      const double* c2;    //            drift      |           amp
      const double* c3;    //            noise_amp  |           (unused)
      const double* amplitude;       // per material, output amplitude of this phase (equilibration or run)
      std::vector<double> white;     // per material, classical open-system thermostat prefactor
      const double* A;               // per material, os-LLG coupling alpha*omega0^4/gamma
      const double* gamma;           // per material, Lorentzian width
      const double* omega0;          // per material, Lorentzian centre
   };
   bool gpu_active();                // the module runs and the GPU takes over its work
   bool gpu_bath();                  // coloured bath for llg-heun on the GPU
   bool gpu_thermostat();            // open-system LLG on the GPU
   void gpu_view(gpu_view_t& view);  // refresh the coefficients for sim::temperature and fill the view

} // end of quantum namespace

#endif //QUANTUM_H_
