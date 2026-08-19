//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SLD_INTERNAL_H_
#define SLD_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the sld module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "sld.hpp"

// sld module headers (pulls in quantum/internal.hpp)
#include "internal.hpp"


namespace sld{

   namespace internal{

      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             set_double_t mass;
             set_double_t V0;
             set_double_t J0;
             set_double_t C0;
             set_double_t damp_lat;
             set_double_t eq_damp_lat;

             set_double_t J0_ms;
             set_double_t C0_ms;
             set_double_t J0_prime;
             set_double_t F_th_sigma;
             set_double_t F_th_sigma_eq;







             // constructor
             mp_t (const unsigned int max_materials = 100){
                mass.set(5.7915e-3);
                V0.set(0.15);
                J0.set(0.904);
                J0_prime.set(3*0.904/7.8);
                J0_ms.set(0.904/2.04028e-23);
                C0.set(0.5);
                C0_ms.set(0.5/2.04028e-23);
                F_th_sigma.set(1.0);
                F_th_sigma_eq.set(1.0);
                damp_lat.set(0.06);
                eq_damp_lat.set(0.6);



             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module
      extern std::vector<sld::internal::mp_t> mp; // array of material properties

      //---------------------------------------------------------------------
      // Quantum noise type selection for both spin and phonon subsystems.
      // Default: sld_classical — backward compatible with classical SLD.
      // Set via input keyword:  spin-lattice:quantum-noise = [classical|quantum|quantum-no-zero]
      //---------------------------------------------------------------------
      enum sld_noise_t {
         sld_classical,           // classical white Gaussian noise (default)
         sld_quantum,             // quantum coth via HO (on-the-fly)
         sld_quantum_no_zero,     // quantum (coth-1) via HO (on-the-fly)
         sld_quantum_fft,         // quantum coth via FFT (pre-generated)
         sld_quantum_no_zero_fft  // quantum (coth-1) via FFT (pre-generated)
      };
      extern sld_noise_t quantum_noise_type;

      //---------------------------------------------------------------------
      // Diagnostic noise export. When enabled, the integrator appends the
      // actual injected noise (spin field added to H_eff and lattice/force
      // noise added to the velocity update) for one atom / one component
      // (x) once per Suzuki-Trotter step to a file. Columns:
      //    time [s]   spin_x   lattice_x
      // Set via:  spin-lattice:export-noise [= filename]
      //           spin-lattice:export-noise-atom = <int>
      //---------------------------------------------------------------------
      extern bool export_noise;                  // enable noise.dat export
      extern std::string export_noise_filename;  // output file (default "noise.dat")
      extern int  export_noise_atom;             // atom index to sample (default 0)
      void write_noise_sample(double t, double spin_x, double lattice_x);

      extern double r_cut_pot; // mechanical potential cutoff
      extern double r_cut_fields; // exchange/coupling cutoff

      extern double dr_init;
      extern double th_velo;

       //for the morse potential
       extern double morse_beta;
       extern double morse_factor;
       extern double alpha_m;
       extern double r0_m;
       extern double morse_D;
       extern bool morse;

      extern bool harmonic; // bool to enable module
      extern bool pseudodipolar;
      extern bool full_neel;


      //extern std::vector<int> sld_neighbour_list_start_index;
      //extern std::vector<int> sld_neighbour_list_end_index;
      //extern std::vector<int> sld_neighbour_list_array;

      extern std::vector<double> x0_coord_array;
      extern std::vector<double> y0_coord_array;
      extern std::vector<double> z0_coord_array;


      extern std::vector <double> x_coord_storage_array;
      extern std::vector <double> y_coord_storage_array;
      extern std::vector <double> z_coord_storage_array;



      extern std::vector<double> forces_array_x;
      extern std::vector<double> forces_array_y;
      extern std::vector<double> forces_array_z;

      extern std::vector<double> fields_array_x;
      extern std::vector<double> fields_array_y;
      extern std::vector<double> fields_array_z;

      extern std::vector<double> velo_array_x;
      extern std::vector<double> velo_array_y;
      extern std::vector<double> velo_array_z;
      extern std::vector<double> potential_eng;
      extern std::vector<double> exch_eng;
      extern std::vector<double> coupl_eng;
      extern std::vector<double> sumJ;
      extern std::vector<double> sumC;

      extern std::vector<int> test_atom_list; //Core atoms of each octant



      void initialise_positions(std::vector<double>& x0_coord_array, // coord vectors for atoms
                  std::vector<double>& y0_coord_array,
                  std::vector<double>& z0_coord_array,
                  std::vector<double>& x_coord_array, // coord vectors for atoms
                  std::vector<double>& y_coord_array,
                  std::vector<double>& z_coord_array,
                  const double dr);

     void thermal_velocity(std::vector<double>& x_velo_array, // coord vectors for atoms
                 std::vector<double>& y_velo_array,
                 std::vector<double>& z_velo_array);

//function to resize vectors and initialise rest of parameters
      void initialise_sld_parameters();

//functions to compute potentials
      void compute_forces_harmonic(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x0_coord_array, // coord vectors for atoms
            const std::vector<double>& y0_coord_array,
            const std::vector<double>& z0_coord_array,
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng);

//functions to compute fields
      void compute_exchange(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            const std::vector<double>& x_spin_array, // coord vectors for atoms
            const std::vector<double>& y_spin_array,
            const std::vector<double>& z_spin_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& fields_array_x, //  vectors for forces
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z);


      void compute_forces_morse(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng);
//
      void compute_sld_coupling(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            const std::vector<double>& x_spin_array, // spin  vectors for atoms
            const std::vector<double>& y_spin_array,
            const std::vector<double>& z_spin_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& fields_array_x, //  vectors for forces
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z);

      void compute_sld_coupling_neel(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            const std::vector<double>& x_spin_array, // spin  vectors for atoms
            const std::vector<double>& y_spin_array,
            const std::vector<double>& z_spin_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& fields_array_x, //  vectors for forces
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z);

//

      void cayley_update(const int start_index,
                  const int end_index,
                  double dt,
                  std::vector<double>& x_spin_array, // coord vectors for atoms
                  std::vector<double>& y_spin_array,
                  std::vector<double>& z_spin_array,
                  std::vector<double>& fields_array_x, //  vectors for fields
                  std::vector<double>& fields_array_y,
                  std::vector<double>& fields_array_z);
   // add spin noise

      void add_spin_noise(const int start_index,
                  const int end_index,
                  double dt,
                  const std::vector<int>& type_array, // type for atom
                  const std::vector<double>& x_spin_array, // coord vectors for atoms
                  const std::vector<double>& y_spin_array,
                  const std::vector<double>& z_spin_array,
                  std::vector<double>& fields_array_x, //  vectors for fields
                  std::vector<double>& fields_array_y,
                  std::vector<double>& fields_array_z,
                  std::vector<double>& Hx_th, //  vectors for fields
                  std::vector<double>& Hy_th,
                  std::vector<double>& Hz_th);

    //MPI variables
    extern std::vector<std::vector<int> > c_octants; //Core atoms of each octant
    extern std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant

    extern std::vector <int> all_atoms_octant_start_index;
    extern std::vector <int> all_atoms_octant_end_index;
    extern std::vector <int> all_atoms_octant;


      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of sld namespace


#endif //SLD_INTERNAL_H_
