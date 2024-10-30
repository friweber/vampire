#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <cstdlib>
#include <random>
#include <complex>
#include <fftw3.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

// Function prototypes
double PSD(const double& omega, const double& T);
std::vector<double> generateNoise(int n, double dt, double T);
void calculate_random_fields(int realizations, int n, double dt, double T);
std::vector<double> noise_array;
void assign_unique_indices(int realizations);




const double A = 152.e5;
const double Gamma = 166.7;
const double omega0 = 282.5;
//const double hbar = 1.05457182E-34;
//const double kB = 1.38064852E-24;
const double hbar = 1;
const double kB = 1;

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void time_series(){

    const int num_atoms = atoms::num_atoms;
    int realizations = num_atoms*3+4;
    std::cout << "Star time series" << std::endl;
    std::cout << realizations << std::endl;

    atoms::x_w_array.resize(atoms::num_atoms, 0.0);
    atoms::y_w_array.resize(atoms::num_atoms, 0.0);
    atoms::z_w_array.resize(atoms::num_atoms, 0.0);
    atoms::x_v_array.resize(atoms::num_atoms, 0.0);
    atoms::y_v_array.resize(atoms::num_atoms, 0.0);
    atoms::z_v_array.resize(atoms::num_atoms, 0.0);

    // check calling of routine if error checking is activated
    if(err::check==true) std::cout << "program::time_series has been called" << std::endl;
    const char* temp_text = std::getenv("SIM_TEMPERATURE");
    double temp = std::atof(temp_text);
    

    std::cout << "Assign directiojs" << std::endl;
    assign_unique_indices(realizations);


    // Resize arrays to the number of atoms

    const int num_steps = static_cast<int>(sim::equilibration_time); // Number of steps
    const double dt = sim::partial_time;

    std::cout << "dt = " << mp::dt << " steps: " << num_steps+1 << std::endl;
    std::cout << "Calculate random fields" << std::endl;
    calculate_random_fields(realizations, num_steps+1, mp::dt, 218.98);
    std::cout << "Calculated random fields" << std::endl;
    atoms::noise_index = 0;

    // //Set initial random!
    // for (int i = 0; i < num_atoms; ++i) {
    //     // Generate random angles for spherical coordinates
    //     double theta = static_cast<double>(std::rand()) / RAND_MAX * 2.0 * M_PI; // 0 to 2*pi
    //     double phi = static_cast<double>(std::rand()) / RAND_MAX * M_PI; // 0 to pi

    //     // Convert spherical coordinates to Cartesian coordinates
    //     atoms::x_spin_array[i] = std::sin(phi) * std::cos(theta);
    //     atoms::y_spin_array[i] = std::sin(phi) * std::sin(theta);
    //     atoms::z_spin_array[i] = std::cos(phi);
    // }

    // Set equilibration temperature only if continue checkpoint not loaded
    if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
    else{
        // Set equilibration temperature
        sim::temperature=218.98;
    }
    // Equilibrate system
    while(sim::time<sim::equilibration_time){
        

    
        sim::integrate(sim::partial_time);
        
        // Calculate magnetisation statistics
        stats::mag_m();

        // Output data
        vout::data();
    }
    // Set temperature and reset stats only if continue checkpoint not loaded
    if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
    else{
        // set simulation temperature
        sim::temperature = temp;

        // Reset mean magnetisation counters
        stats::mag_m_reset();
    }


    // // Vector to hold temperature data from file
    // std::vector<double> temperature_data;
    // std::ifstream file(std::getenv("TEMPERATURE_FILE")); // File containing temperature data
    // double temp_list;

    // // Read the temperature data from file
    // while (file >> temp_list) {
    //     temperature_data.push_back(temp_list);
    // }
    // file.close();

    // // Ensure there is temperature data
    // if (temperature_data.empty()) {
    //     std::cerr << "Error: semperature data loaded." << std::endl;
    //     exit(EXIT_FAILURE); // Exit the program if no data is loaded
    // }

    std::cout << "dt = " << mp::dt << " steps: " << num_steps+1 << std::endl;
    std::cout << "Calculate random fields" << std::endl;
    const int num_steps_2 = static_cast<int>(sim::total_time); 
    calculate_random_fields(realizations, num_steps_2+1, mp::dt, temp);
    std::cout << "Calculated random fields" << std::endl;

    atoms::noise_index = 0;
    // Perform Time Series
    while (sim::time < sim::equilibration_time + sim::total_time) {

        // Integrate system
        sim::integrate(sim::partial_time);
        
        // Calculate magnetisation statistics
        stats::mag_m();

        // Output data
        vout::data();
    }
}


void assign_unique_indices(int realizations) {
    const int num_atoms = atoms::num_atoms;
    const int indices_per_atom = 3; // x, y, z

    // Ensure we have enough unique indices for all atoms
    if (realizations < num_atoms * indices_per_atom) {
        throw std::runtime_error("Not enough realizations for unique indices");
    }

    std::vector<int> all_indices(realizations);
    std::iota(all_indices.begin(), all_indices.end(), 0); // Fill with 0, 1, ..., realizations-1

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(all_indices.begin(), all_indices.end(), gen);

    atoms::atom_idx_x.resize(num_atoms);
    atoms::atom_idx_y.resize(num_atoms);
    atoms::atom_idx_z.resize(num_atoms);

    for (int atom = 0; atom < num_atoms; atom++) {
        int base_index = atom * indices_per_atom;
        atoms::atom_idx_x[atom] = all_indices[base_index];
        atoms::atom_idx_y[atom] = all_indices[base_index + 1];
        atoms::atom_idx_z[atom] = all_indices[base_index + 2];
    }
}

double PSD(const double& omega, const double& T) {
    if (omega == 0.0)
        return 0.0;

    double x = hbar * omega / (2 * kB * T);  // hbar and kB constants
    double lorentzian = A * Gamma * hbar * kB * omega / ((omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega);
    double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);  // Stabilize coth calculation near zero
    return (coth - 1) * lorentzian; //No Zero Quantum Noise
    //return (coth - 1); //Non lorenzian no zero Quantum Noise
    //return coth * lorentzian; //Quantum Noies
    //return 1.0/x * lorentzian; //Classical Noise
}


std::vector<double> generateNoise(int N, double dt, double T) {
    // Allocate memory
    double *in = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));
    double *result = (double*) fftw_malloc(sizeof(double) * N);
    std::vector<double> noise(N);

    // Generate white noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt));
    for (int i = 0; i < N; ++i) {
        in[i] = dist(gen);
    }

    // Create FFT plan and execute
    fftw_plan forward = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);

    // Calculate frequency step
    double df = 1.0 / (N * dt);

    // Apply PSD
    for (int i = 0; i <= N/2; i++) {
        double freq = i * df;
        double omega = freq * 2 * M_PI;
        double magnitude = std::sqrt(PSD(omega, T));
        out[i][0] *= magnitude;
        out[i][1] *= magnitude;
    }

    // Inverse FFT
    fftw_plan backward = fftw_plan_dft_c2r_1d(N, out, result, FFTW_ESTIMATE);
    fftw_execute(backward);

    // Normalize
    for (int i = 0; i < N; i++) {
        result[i] /= N;
    }

    for (int i = 0; i < N; i++) {
        noise[i] = result[i];
    }

    // Clean up
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(in);
    fftw_free(out);
    fftw_free(result);

    return noise;
}


void calculate_random_fields(int realizations, int n, double dt, double T) {
    // Resize atoms::noise_field to hold noise for each atom
    atoms::noise_field.resize(realizations, std::vector<double>(n));

    const int bar_width = 50;
    int last_printed_percent = -1;
    
    for (int r = 0; r < realizations; ++r) {
        std::vector<double> noise = generateNoise(n, dt, T);
        for (int j = 0; j < n; ++j) {
            atoms::noise_field[r][j] = noise[j];
        }

        // Calculate current progress percentage
        int current_percent = static_cast<int>(std::round((r + 1) * 100.0 / realizations));

        // Update loading bar if progress has increased by at least 1%
        if (current_percent > last_printed_percent) {
            float progress = static_cast<float>(r + 1) / realizations;
            int pos = static_cast<int>(bar_width * progress);

            std::cout << "\r[";
            for (int i = 0; i < bar_width; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << current_percent << "%";
            std::cout.flush();

            last_printed_percent = current_percent;
        }
    }
    std::cout << std::endl;  // Move to the next line after completion
}

} // end of namespace program
