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
//   Public interface for the quantum thermostat module.
//
//------------------------------------------------------------------------------

#ifndef QUANTUM_H_
#define QUANTUM_H_

// C++ standard library headers
#include <cstdint>
#include <string>
#include <vector>

namespace quantum{

   //------------------------------------------------------------------------
   // Module initialization
   //------------------------------------------------------------------------
   void initialize();

   //------------------------------------------------------------------------
   // Release persistent resources (FFTW plans + buffers). Called once at
   // program shutdown. Safe to call even if the module was never enabled.
   //------------------------------------------------------------------------
   void cleanup();

   //------------------------------------------------------------------------
   // Integrate quantum LLG equation (one time step)
   //------------------------------------------------------------------------
   void llg();

   //------------------------------------------------------------------------
   // Get FFT noise field for a given atom and component (0=x, 1=y, 2=z) at
   // the current fine-step index. One sample per RK4 step — the same value
   // is reused across all four sub-stages, matching the LSF_RK4 / HO convention.
   //------------------------------------------------------------------------
   double get_field(int atom, int component);

   //------------------------------------------------------------------------
   // Advance fine-grained noise time index by one coarse step
   //------------------------------------------------------------------------
   void increment_time();

   //------------------------------------------------------------------------
   // Input file parameter matching
   //------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word,
                              std::string const value, std::string const unit,
                              int const line);

   //------------------------------------------------------------------------
   // Material file parameter matching
   //------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value,
                                 std::string const unit, int const line,
                                 int const super_index, const int sub_index);

   //------------------------------------------------------------------------
   // Public SLD-noise API
   //
   // The spin-lattice (sld) module uses the quantum module's coloured-noise
   // bath through this small surface only — it never touches quantum::internal.
   // Amplitudes are derived entirely from SLD material parameters passed in
   // via initialize(), so an SLD quantum-noise run needs no quantum-lorentzian-*
   // inputs (only the universal coth / coth-1 spectral shape is reused).
   //------------------------------------------------------------------------
   namespace sld_noise{

      /// Noise kind selected by spin-lattice:noise-type
      enum class kind_t {
         classical,
         quantum,              // coth spectral shape via HO (on-the-fly)
         quantum_no_zero,      // (coth-1) spectral shape via HO (on-the-fly)
         quantum_fft,          // coth spectral shape via FFT (pre-generated)
         quantum_no_zero_fft   // (coth-1) spectral shape via FFT (pre-generated)
      };

      /// Per-material SLD parameters needed to derive the noise amplitudes,
      /// supplied by the sld module (keeps quantum free of sld internals).
      struct material_params {
         double mass;        // atomic mass [kg] (0 = ASD-only, skip phonon bath)
         double damp_lat;    // lattice (phonon) damping constant
         double V0;          // harmonic spring constant [eV/Ang^2]
         double H_th_sigma;  // spin thermal-field prefactor (sigma, before *sqrt(T))
      };

      /// One-time setup. classical -> no-op (white noise stays in sld).
      /// For FFT variants, n_fine is the total number of integration steps
      /// (equilibration + production); used to size the pre-generated array.
      void initialize(kind_t kind, int num_atoms,
                      const std::vector<material_params>& mats,
                      uint64_t n_fine = 0);

      /// Draw one spin + one phonon sample per atom for the current step.
      /// For FFT kinds, advances the internal step index; no FFT work done here.
      void generate(int num_atoms);

      /// Per-atom noise accessors (component 0=x, 1=y, 2=z), valid after generate().
      double spin(int atom, int component);    // added to H_eff
      double phonon(int atom, int component);  // added to the force/velocity update

      /// Enable per-step export of the injected spin noise to a text file.
      /// Each call to export_spin_noise_step() appends one line:
      ///   time[s] \t spin_x \t spin_y \t spin_z
      /// for the chosen atom. Call once after initialize() when export is wanted.
      void enable_spin_noise_export(const std::string& filename, int atom);

      /// Append one row to the export file.  No-op if export not enabled.
      void export_spin_noise_step(double t);

   } // end of sld_noise namespace

} // end of quantum namespace

#endif // QUANTUM_H_
