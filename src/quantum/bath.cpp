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
#include <cmath>
#include <iostream>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      //---------------------------------------------------------------------------
      // Function to set up one bath: store its parameters, allocate the per-site
      // state and build the coefficients for the current temperature. No random
      // numbers are drawn here; the generator is not seeded until sim::run().
      //---------------------------------------------------------------------------
      void setup_bath(bath_t& b, const spectrum_t spectrum, const generator_t generator,
                      const int n_sites, const double dt, const double T_scale,
                      const std::vector<double>& amp, const std::vector<double>& amp_eq,
                      const bool lorentzian, const uint64_t n_run){

         b.spectrum = spectrum;
         b.generator = generator;
         b.lorentzian = lorentzian;
         b.n_sites = n_sites;
         b.n_modes = n_bath_modes;
         b.dt = dt;
         b.T_scale = T_scale;
         b.amp = amp;
         b.amp_eq = amp_eq;

         b.x.assign(n_sites, 0.0);
         b.y.assign(n_sites, 0.0);
         b.z.assign(n_sites, 0.0);

         // pre-generated noise: shape the window filter, allocate the buffers
         if(generator == pre_generated){
            fft_setup(b, n_run);
            return;
         }

         // white noise needs no state (thermostat only)
         if(spectrum == classical) return;

         // auxiliary modes, one set per site and component
         const size_t n_state = static_cast<size_t>(n_sites)*b.n_modes;
         b.s_x.assign(n_state, 0.0);
         b.s_y.assign(n_state, 0.0);
         b.s_z.assign(n_state, 0.0);
         b.scratch.assign(3*b.n_modes + 3, 0.0);

         const double T_scaled = sim::temperature*b.T_scale;

         if(spectrum == quantum_zero){
            ou_coefficients(b, T_scaled);
         }
         else{
            log_bath_fit(b, T_scaled);
            log_bath_coefficients(b, T_scaled);
         }

         b.last_T_scaled = T_scaled;

         return;

      }

      //---------------------------------------------------------------------------
      // Function to rebuild the temperature dependent coefficients when
      // sim::temperature has moved (on-the-fly generator only)
      //---------------------------------------------------------------------------
      void refresh_temperature(bath_t& b){

         if(b.generator != on_the_fly || b.spectrum == classical) return;

         const double T_scaled = sim::temperature*b.T_scale;
         if(T_scaled == b.last_T_scaled) return;

         if(b.spectrum == quantum_zero) ou_coefficients(b, T_scaled);
         else log_bath_coefficients(b, T_scaled);

         b.last_T_scaled = T_scaled;

         return;

      }

      //---------------------------------------------------------------------------
      // Function to draw this step's samples for every site of a bath
      //---------------------------------------------------------------------------
      void draw(bath_t& b){

         // equilibration uses its own amplitudes, as the integrators do for damping
         const std::vector<double>& amplitude = (sim::time < sim::equilibration_time) ? b.amp_eq : b.amp;

         if(b.generator == pre_generated){
            fft_draw(b, amplitude);
            return;
         }

         if(b.spectrum == classical){
            draw_white_thermostat(b);
            return;
         }

         refresh_temperature(b);

         if(b.spectrum == quantum_zero){
            for(int site = 0; site < b.n_sites; site++){
               ou_draw(b, site, amplitude[atoms::type_array[site]]);
            }
         }
         else{
            for(int site = 0; site < b.n_sites; site++){
               log_bath_draw(b, site, amplitude[atoms::type_array[site]]);
            }
         }

         return;

      }

      //---------------------------------------------------------------------------
      // Function to draw white noise for the thermostat. The prefactor
      // sqrt(2 Gamma A T/(S0 dt)) balances the Lorentzian damping, and the
      // temperature carries the per-material rescaling that the classical
      // thermal field applies.
      //---------------------------------------------------------------------------
      void draw_white_thermostat(bath_t& b){

         std::vector<double> T_scaled_mat(mp::material.size());
         for(size_t m = 0; m < mp::material.size(); m++){
            double T = sim::temperature;
            if(sim::local_temperature) T = mp::material[m].temperature;
            const double a  = mp::material[m].temperature_rescaling_alpha;
            const double Tc = mp::material[m].temperature_rescaling_Tc;
            const double T_rescaled = (T < Tc) ? Tc*std::pow(T/Tc, a) : T;
            T_scaled_mat[m] = T_rescaled*b.T_scale;
         }

         for(int atom = 0; atom < b.n_sites; atom++){
            const int m = atoms::type_array[atom];
            const double prefactor = std::sqrt(2.0*material_gamma[m]*material_A[m]*T_scaled_mat[m]/(material_S0[m]*b.dt));
            b.x[atom] = prefactor*mtrandom::gaussian();
            b.y[atom] = prefactor*mtrandom::gaussian();
            b.z[atom] = prefactor*mtrandom::gaussian();
         }

         return;

      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Function to report whether a coloured bath replaces the thermal field
   //---------------------------------------------------------------------------
   bool enabled(){
      return internal::mode == internal::bath;
   }

   //---------------------------------------------------------------------------
   // Function to draw this step's samples for the spin and lattice baths
   //---------------------------------------------------------------------------
   void generate(){

      using namespace internal;

      if(mode != bath) return;

      draw(spin_bath);
      if(lattice_bath.n_sites > 0) draw(lattice_bath);

      if(export_noise) export_sample(spin_bath.x[export_atom], spin_bath.y[export_atom], spin_bath.z[export_atom]);

      return;

   }

   //---------------------------------------------------------------------------
   // Function to return this step's spin-bath sample for one atom (Tesla)
   //---------------------------------------------------------------------------
   double field(const int atom, const int component){
      if(component == 0) return internal::spin_bath.x[atom];
      if(component == 1) return internal::spin_bath.y[atom];
      return internal::spin_bath.z[atom];
   }

   //---------------------------------------------------------------------------
   // Function to add the spin-bath sample to the external field of a range of atoms
   //---------------------------------------------------------------------------
   void add_field(const int start_index, const int end_index){

      using namespace internal;

      if(mode != bath) return;

      for(int atom = start_index; atom < end_index; atom++){
         atoms::x_total_external_field_array[atom] += spin_bath.x[atom];
         atoms::y_total_external_field_array[atom] += spin_bath.y[atom];
         atoms::z_total_external_field_array[atom] += spin_bath.z[atom];
      }

      return;

   }

   //---------------------------------------------------------------------------
   // Function to store the lattice parameters handed over by the spin-lattice
   // module. Called from sld::initialize(), before quantum::initialize().
   //---------------------------------------------------------------------------
   void set_lattice_parameters(const std::vector<double>& mass,
                               const std::vector<double>& damping,
                               const std::vector<double>& damping_eq){

      internal::lattice_mass = mass;
      internal::lattice_damping = damping;
      internal::lattice_damping_eq = damping_eq;
      internal::lattice_parameters_set = true;

      return;

   }

   //---------------------------------------------------------------------------
   // Function to return this step's lattice-bath sample for one atom (force)
   //---------------------------------------------------------------------------
   double lattice_field(const int atom, const int component){
      if(component == 0) return internal::lattice_bath.x[atom];
      if(component == 1) return internal::lattice_bath.y[atom];
      return internal::lattice_bath.z[atom];
   }

} // end of quantum namespace
