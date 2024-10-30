//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Curie Temperature program
///
/// @details Performs a temperature loop to calculate temperature dependent magnetisation
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	09/03/2011
///=====================================================================================
///

// Standard Libraries
#include <iostream>

#include <cstdlib>
#include <random>
#include <complex>
#include <fftw3.h>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"



namespace program{

// Function prototypes
double PSD(double omega, double T);
void fft(fftw_complex* in, fftw_complex* out, int n);
void ifft(fftw_complex* in, fftw_complex* out, int n);
std::vector<std::complex<double>> generateNoise(int n, double dt, double T);
void calculate_random_fields(int realizations, int n, double dt, double T);
void assign_unique_indices(int realizations);
const int realizations = 500;

const  double A = 152.e5;
const double Gamma = 166.7;
const double omega0 = 282.5;
//const double hbar = 1.05457182E-34;
//const double kB = 1.38064852E-24;
const double hbar = 1;
const double kB = 1;


/// @brief Function to calculate the temperature dependence of the magnetisation
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised
/// accoring to the input flag - either randomly or ordered.For the ordered case the temperature sequence
/// increases from zero, for the random case the temperature decreases from the maximum temperature. After
/// initialisation the sytem is equilibrated for sim::equilibration timesteps.
///
/// @section notes Implementation Notes
/// Capable of hot>cold or cold>hot calculation.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
///
/// @param[in] init Determines whether the system is initialised randomly (0) or ordered (1)
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		11/01/2010
///	Revision:	09/03/2011
///=====================================================================================
///

int curie_temperature(){
	std::cout << "Resize" << std::endl;
	atoms::x_w_array.resize(atoms::num_atoms, 0.0);
    atoms::y_w_array.resize(atoms::num_atoms, 0.0);
    atoms::z_w_array.resize(atoms::num_atoms, 0.0);
    atoms::x_v_array.resize(atoms::num_atoms, 0.0);
    atoms::y_v_array.resize(atoms::num_atoms, 0.0);
    atoms::z_v_array.resize(atoms::num_atoms, 0.0);

	std::cout << "START CURIE" << std::endl;
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::curie_temperature has been called" << std::endl;}

	// Set starting temperature
   // Initialise sim::temperature
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){
      sim::temperature+=sim::delta_temperature;
   }
   else sim::temperature=sim::Tmin;


	const int num_atoms = atoms::num_atoms;
    // Resize arrays to the number of atoms
    atoms::random_phi_array.resize(num_atoms+2);
    atoms::random_theta_array.resize(num_atoms+2);

    // Random number generator setup
    std::mt19937 gen(12345);  // Random number generator with a fixed seed
    std::uniform_real_distribution<> dist_p(0.0, 2 * M_PI);
    std::uniform_real_distribution<> dist_t(0.0, M_PI);
    // Assign random rotations to each atom
    for (int atom = 0; atom < num_atoms; atom++) {
        atoms::random_phi_array[atom] = dist_p(gen);  
        atoms::random_theta_array[atom] = dist_t(gen);
    }
    std::cout << "Assign directions" << std::endl;
    assign_unique_indices(realizations);

    const int num_steps = static_cast<int>(sim::equilibration_time+sim::loop_time); // Number of steps
    const double dt = sim::partial_time;

    std::cout << "dt = " << mp::dt << " steps: " << num_steps+1 << std::endl;
    std::cout << "Calculate random fields" << std::endl;
    calculate_random_fields(realizations, num_steps+1, mp::dt, sim::temperature);


	// Perform Temperature Loop
	while(sim::temperature<=sim::Tmax){

		// Equilibrate system
		sim::integrate(sim::partial_time);

		// Reset mean magnetisation counters
		stats::mag_m_reset();

		// Reset start time
		int start_time=sim::time;

		// Simulate system
		while(sim::time<sim::loop_time+start_time){

			// Integrate system
			sim::integrate(sim::partial_time);

			// Calculate magnetisation statistics
			stats::mag_m();

		}

		// Output data
		vout::data();

		// Increment temperature
		sim::temperature+=sim::delta_temperature;

	} // End of temperature loop

	return EXIT_SUCCESS;
}
// void assign_unique_indices(int realizations) {
//     const int num_atoms = atoms::num_atoms;
//     std::vector<int> indices(realizations);
//     std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, ..., realizations-1
//     atoms::atom_idx_x.resize(num_atoms);
//     atoms::atom_idx_y.resize(num_atoms);
//     atoms::atom_idx_z.resize(num_atoms);
//     std::random_device rd;
//     std::mt19937 gen(rd());

//     for (int atom = 0; atom < num_atoms; atom++) {
//         std::shuffle(indices.begin(), indices.end(), gen);

//         // Assign first three unique indices to atom
//         atoms::atom_idx_x[atom] = indices[0];
//         atoms::atom_idx_y[atom] = indices[1];
//         atoms::atom_idx_z[atom] = indices[2];
//     }
// }

// double PSD(double omega, double T) {
//     if (omega == 0.0)
//         return 0.0;

//     double x = hbar * omega / (2 * kB * T);
//     double lorentzian = A * Gamma * hbar * kB * omega / ((omega0 * omega0 - omega * omega) * (omega0 * omega * omega0 - omega * omega) + Gamma * Gamma * omega * omega);
//     double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);  // Stabilize coth calculation near zero
//     return (coth - 1) * lorentzian;
// }

// void fft(fftw_complex* in, fftw_complex* out, int n) {
//     fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
// }

// void ifft(fftw_complex* in, fftw_complex* out, int n) {
//     fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (int i = 0; i < n; i++) {
//         out[i][0] /= n;
//         out[i][1] /= n;
//     }
// }

// std::vector<std::complex<double>> generateNoise(int n, double dt, double T) {
//     fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
//     fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
//     fftw_complex *convolution = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
//     fftw_complex *final_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

//     if (!in || !out || !convolution || !final_result) {
//         if (in) fftw_free(in);
//         if (out) fftw_free(out);
//         if (convolution) fftw_free(convolution);
//         if (final_result) fftw_free(final_result);
//         exit(EXIT_FAILURE);
//     }

//     std::vector<std::complex<double>> noise(n);
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt));

//     for (int i = 0; i < n; ++i) {
//         in[i][0] = dist(gen);
//         in[i][1] = dist(gen);
//     }

//     // Perform FFT on the input function
//     fft(in, out, n);

//     // Calculate frequency step
//     double df = 1.0 / (n * dt);

//     // Symmetrize PSD and convolution
//     for (int i = 0; i < n / 2; i++) {
//         double freq = i * df;
//         double magnitude = PSD(2 * M_PI * freq, T);  // Take sqrt for amplitude spectral density

//         // Check for NaN magnitude and replace if necessary
//         if (std::isnan(magnitude)) {
//             magnitude = 0.0; // Replace NaN with zero or a small value
//         }

//         convolution[i][0] = out[i][0] * sqrt(magnitude);
//         convolution[i][1] = out[i][1] * sqrt(magnitude);

//         // Symmetrize the values for negative frequencies
//         int neg_index = n - i - 1;
//         convolution[neg_index][0] = convolution[i][0];
//         convolution[neg_index][1] = convolution[i][1];

//         // Check for NaN in the convolution result and replace if necessary
//         if (std::isnan(convolution[i][0])) {
//             convolution[i][0] = 0.0; // Replace NaN with zero or a small value
//             convolution[neg_index][0] = 0.0; // Ensure symmetry
//         }
//         if (std::isnan(convolution[i][1])) {
//             convolution[i][1] = 0.0; // Replace NaN with zero or a small value
//             convolution[neg_index][1] = 0.0; // Ensure symmetry
//         }
//     }

//     // Perform inverse FFT on the result
//     ifft(convolution, final_result, n);

//     for (int i = 0; i < n; ++i) {
//         noise[i] = std::complex<double>(final_result[i][0], final_result[i][1]);
//     }

//     // Free memory
//     fftw_free(in);
//     fftw_free(out);
//     fftw_free(convolution);
//     fftw_free(final_result);

//     return noise;
// }


// void calculate_random_fields(int realizations, int n, double dt, double T) {
//     // Resize atoms::noise_field to hold noise for each atom
//     atoms::noise_field.resize(realizations, std::vector<double>(n));
//     for (int r = 0; r < realizations; ++r) {
//         std::vector<std::complex<double>> noise = generateNoise(n, dt, T);
//         for (int j = 0; j < n; ++j) {
//             int idx_x = atoms::atom_idx_x[0];
//             if (r == idx_x){
// 			    //std::cout << noise[j].real() << std::endl;
// 		    }
//             atoms::noise_field[r][j] = noise[j].real();
//         }
//     }
// }

} // end of namespace program
