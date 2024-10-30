/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
/// @internal
/// Created: 05/02/2011
/// Revision: ---
///=====================================================================================
///

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <functional>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"

//Function prototypes
int calculate_spin_fields(const int, const int);
int calculate_external_fields(const int, const int);

namespace LLG_arrays {
    // Local arrays for LLG integration
    std::vector <double> x_heun_array;
	std::vector <double> y_heun_array;
	std::vector <double> z_heun_array;

    std::vector<double> x_euler_array;
    std::vector<double> y_euler_array;
    std::vector<double> z_euler_array;

    std::vector<double> x_spin_storage_array;
    std::vector<double> y_spin_storage_array;
    std::vector<double> z_spin_storage_array;

    std::vector<double> x_initial_spin_array;
    std::vector<double> y_initial_spin_array;
    std::vector<double> z_initial_spin_array;

    std::vector<double> w_x_initial_spin_array;
    std::vector<double> w_y_initial_spin_array;
    std::vector<double> w_z_initial_spin_array;

    std::vector<double> v_x_initial_spin_array;
    std::vector<double> v_y_initial_spin_array;
    std::vector<double> v_z_initial_spin_array;

    std::vector<double> k1_x_array;
    std::vector<double> k1_y_array;
    std::vector<double> k1_z_array;
    std::vector<double> k1_vx_array;
    std::vector<double> k1_vy_array;
    std::vector<double> k1_vz_array;
    std::vector<double> k1_wx_array;
    std::vector<double> k1_wy_array;
    std::vector<double> k1_wz_array;

    std::vector<double> k2_x_array;
    std::vector<double> k2_y_array;
    std::vector<double> k2_z_array;
    std::vector<double> k2_vx_array;
    std::vector<double> k2_vy_array;
    std::vector<double> k2_vz_array;
    std::vector<double> k2_wx_array;
    std::vector<double> k2_wy_array;
    std::vector<double> k2_wz_array;

    std::vector<double> k3_x_array;
    std::vector<double> k3_y_array;
    std::vector<double> k3_z_array;
    std::vector<double> k3_vx_array;
    std::vector<double> k3_vy_array;
    std::vector<double> k3_vz_array;
    std::vector<double> k3_wx_array;
    std::vector<double> k3_wy_array;
    std::vector<double> k3_wz_array;

    std::vector<double> k4_x_array;
    std::vector<double> k4_y_array;
    std::vector<double> k4_z_array;
    std::vector<double> k4_vx_array;
    std::vector<double> k4_vy_array;
    std::vector<double> k4_vz_array;
    std::vector<double> k4_wx_array;
    std::vector<double> k4_wy_array;
    std::vector<double> k4_wz_array;

    bool LLG_set = false; ///< Flag to define state of LLG arrays (initialized/uninitialized)
}

namespace sim {
    int LLGinit() {
        // Check calling of routine if error checking is activated
        if (err::check == true) { std::cout << "sim::LLG_init has been called" << std::endl; }

        using namespace LLG_arrays;

        std::cout << "Non Markovian LLG has been called" << std::endl;

        // Initial spin arrays for components of s
        x_initial_spin_array.resize(atoms::num_atoms, 0.0);
        y_initial_spin_array.resize(atoms::num_atoms, 0.0);
        z_initial_spin_array.resize(atoms::num_atoms, 0.0);

        // Initial spin arrays for components of w
        w_x_initial_spin_array.resize(atoms::num_atoms, 0.0);
        w_y_initial_spin_array.resize(atoms::num_atoms, 0.0);
        w_z_initial_spin_array.resize(atoms::num_atoms, 0.0);

        // Initial spin arrays for components of v
        v_x_initial_spin_array.resize(atoms::num_atoms, 0.0);
        v_y_initial_spin_array.resize(atoms::num_atoms, 0.0);
        v_z_initial_spin_array.resize(atoms::num_atoms, 0.0);

        // Allocate memory for RK4 intermediate steps
        k1_x_array.resize(atoms::num_atoms, 0.0);
        k1_y_array.resize(atoms::num_atoms, 0.0);
        k1_z_array.resize(atoms::num_atoms, 0.0);
        k1_vx_array.resize(atoms::num_atoms, 0.0);
        k1_vy_array.resize(atoms::num_atoms, 0.0);
        k1_vz_array.resize(atoms::num_atoms, 0.0);
        k1_wx_array.resize(atoms::num_atoms, 0.0);
        k1_wy_array.resize(atoms::num_atoms, 0.0);
        k1_wz_array.resize(atoms::num_atoms, 0.0);

        k2_x_array.resize(atoms::num_atoms, 0.0);
        k2_y_array.resize(atoms::num_atoms, 0.0);
        k2_z_array.resize(atoms::num_atoms, 0.0);
        k2_vx_array.resize(atoms::num_atoms, 0.0);
        k2_vy_array.resize(atoms::num_atoms, 0.0);
        k2_vz_array.resize(atoms::num_atoms, 0.0);
        k2_wx_array.resize(atoms::num_atoms, 0.0);
        k2_wy_array.resize(atoms::num_atoms, 0.0);
        k2_wz_array.resize(atoms::num_atoms, 0.0);

        k3_x_array.resize(atoms::num_atoms, 0.0);
        k3_y_array.resize(atoms::num_atoms, 0.0);
        k3_z_array.resize(atoms::num_atoms, 0.0);
        k3_vx_array.resize(atoms::num_atoms, 0.0);
        k3_vy_array.resize(atoms::num_atoms, 0.0);
        k3_vz_array.resize(atoms::num_atoms, 0.0);
        k3_wx_array.resize(atoms::num_atoms, 0.0);
        k3_wy_array.resize(atoms::num_atoms, 0.0);
        k3_wz_array.resize(atoms::num_atoms, 0.0);

        k4_x_array.resize(atoms::num_atoms, 0.0);
        k4_y_array.resize(atoms::num_atoms, 0.0);
        k4_z_array.resize(atoms::num_atoms, 0.0);
        k4_vx_array.resize(atoms::num_atoms, 0.0);
        k4_vy_array.resize(atoms::num_atoms, 0.0);
        k4_vz_array.resize(atoms::num_atoms, 0.0);
        k4_wx_array.resize(atoms::num_atoms, 0.0);
        k4_wy_array.resize(atoms::num_atoms, 0.0);
        k4_wz_array.resize(atoms::num_atoms, 0.0);

        LLG_set = true;

        return EXIT_SUCCESS;
    }

    // Define cross product of two 3D vectors
    std::vector<double> cross(const std::vector<double>& a, const std::vector<double>& b) {
        return {
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        };
    }

    std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b) {
        return {
            a[0]+b[0],
            a[1]+b[1],
            a[2]+b[2],
        };
    }

    std::vector<double> scalar(const double& a, const std::vector<double>& b) {
        return {
            a*b[0],
            a*b[1],
            a*b[2],
        };
    }

    // ODE system for spin dynamics
    std::vector<double> spinDynamics(const std::vector<double>& y, const std::vector<double>& b) {
        const double A = 152.e5;
        const double Gamma = 166.7;
        const double omega0 = 282.5;
        const double gyromagnetic_ratio = 1.;

        // Unpack y into S, V, W
        std::vector<double> S = {y[0], y[1], y[2]};
        std::vector<double> V = {y[3], y[4], y[5]};
        std::vector<double> W = {y[6], y[7], y[8]};

        // Compute dS/dt
        std::vector<double> dSdt = {
            gyromagnetic_ratio * (S[1]*(b[2]+V[2]) - S[2]*(b[1]+V[1])),
            gyromagnetic_ratio * (S[2]*(b[0]+V[0]) - S[0]*(b[2]+V[2])), 
            gyromagnetic_ratio * (S[0]*(b[1]+V[1]) - S[1]*(b[0]+V[0]))
        };

        // Compute dV/dt and dW/dt
        std::vector<double> dVdt = W;
        std::vector<double> dWdt = {
            -omega0*omega0 * V[0] - Gamma * W[0] + A * gyromagnetic_ratio * S[0],
            -omega0*omega0 * V[1] - Gamma * W[1] + A * gyromagnetic_ratio * S[1],
            -omega0*omega0 * V[2] - Gamma * W[2] + A * gyromagnetic_ratio * S[2]
        };

        // Return the concatenated derivatives
        std::vector<double> dydt = {dSdt[0], dSdt[1], dSdt[2], dVdt[0], dVdt[1], dVdt[2], dWdt[0], dWdt[1], dWdt[2]};
        return dydt;
    }

    int LLG_Heun() {
        // Check calling of routine if error checking is activated
        if (err::check == true) { std::cout << "sim::LLG_RK4 has been called" << std::endl; }
        using namespace LLG_arrays;

        // Check for initialization of LLG integration arrays
        if (LLG_set == false) sim::LLGinit();

        // Local variables for system integration
        const int num_atoms = atoms::num_atoms;
        std::vector<double> y(9);
        std::vector<double> y_temp(9);
        double S_new[3];    // New Local Spin Moment
        double V_new[3];    // New Local V Moment
        double W_new[3];    // New Local W Moment
        double mod_S;       // Magnitude of spin moment

        
        // Store initial spin positions
        for (int atom = 0; atom < num_atoms; atom++) {
            x_initial_spin_array[atom] = atoms::x_spin_array[atom];
            y_initial_spin_array[atom] = atoms::y_spin_array[atom];
            z_initial_spin_array[atom] = atoms::z_spin_array[atom];

            w_x_initial_spin_array[atom] = atoms::x_w_array[atom];
            w_y_initial_spin_array[atom] = atoms::y_w_array[atom];
            w_z_initial_spin_array[atom] = atoms::z_w_array[atom];

            v_x_initial_spin_array[atom] = atoms::x_v_array[atom];
            v_y_initial_spin_array[atom] = atoms::y_v_array[atom];
            v_z_initial_spin_array[atom] = atoms::z_v_array[atom];
        }

    

        // Calculate external fields only once before RK4 steps
        calculate_external_fields(0, num_atoms);
        calculate_spin_fields(0, num_atoms);


        // const double S_x = atoms::x_spin_array[0];
        // const double S_y = atoms::y_spin_array[0];
        // const double S_z = atoms::z_spin_array[0];

        // const double V_x = atoms::x_v_array[0];
        // const double V_y = atoms::y_v_array[0];
        // const double V_z = atoms::z_v_array[0];

        // const double W_x = atoms::x_w_array[0];
        // const double W_y = atoms::y_w_array[0];
        // const double W_z = atoms::z_w_array[0];



        

        for (int atom = 0; atom < num_atoms; ++atom) {
            // Store local field in H
            std::vector<double> H = {
                atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
            };
            // Prepare initial state y
            y[0] = atoms::x_spin_array[atom];
            y[1] = atoms::y_spin_array[atom];
            y[2] = atoms::z_spin_array[atom];
            y[3] = atoms::x_v_array[atom];
            y[4] = atoms::y_v_array[atom];
            y[5] = atoms::z_v_array[atom];
            y[6] = atoms::x_w_array[atom];
            y[7] = atoms::y_w_array[atom];
            y[8] = atoms::z_w_array[atom];

            

    

            

            // Calculate k1
            auto k1 = spinDynamics(y, H);
            k1_x_array[atom] = k1[0];
            k1_y_array[atom] = k1[1];
            k1_z_array[atom] = k1[2];
            k1_vx_array[atom] = k1[3];
            k1_vy_array[atom] = k1[4];
            k1_vz_array[atom] = k1[5];
            k1_wx_array[atom] = k1[6];
            k1_wy_array[atom] = k1[7];
            k1_wz_array[atom] = k1[8];

            // Update temporary state for k2
            for (size_t i = 0; i < y.size(); ++i) {
                y_temp[i] = y[i] + 0.5 * mp::dt * k1[i];
            }

            // Update arrays with temporary values
            atoms::x_spin_array[atom] = y_temp[0];
            atoms::y_spin_array[atom] = y_temp[1];
            atoms::z_spin_array[atom] = y_temp[2];
            atoms::x_v_array[atom] = y_temp[3];
            atoms::y_v_array[atom] = y_temp[4];
            atoms::z_v_array[atom] = y_temp[5];
            atoms::x_w_array[atom] = y_temp[6];
            atoms::y_w_array[atom] = y_temp[7];
            atoms::z_w_array[atom] = y_temp[8];
        }

        // Recalculate spin fields after k1 step
        calculate_spin_fields(0, num_atoms);



        // const double H_ext_x = atoms::x_total_external_field_array[0];
        // const double H_ext_y = atoms::y_total_external_field_array[0];
        // const double H_ext_z = atoms::z_total_external_field_array[0];

        // const double H_spin_x = atoms::x_total_spin_field_array[0];
        // const double H_spin_y = atoms::y_total_spin_field_array[0];
        // const double H_spin_z = atoms::z_total_spin_field_array[0];


        // std::cout << std::setprecision(17);
        // std::cout << S_x << " " << S_y << " " << S_z << " ";
        // std::cout << V_x << " " << V_y << " " << V_z << " ";
        // std::cout << W_x << " " << W_y << " " << W_z << " ";
        // std::cout << H_ext_x << " " << H_ext_y << " " << H_ext_z << " ";
        // std::cout << H_spin_x << " " << H_spin_y << " " << H_spin_z << std::endl;

        // Second RK4 Step (k2)
        for (int atom = 0; atom < num_atoms; ++atom) {
            // Prepare state y_temp from previous step
            y_temp[0] = atoms::x_spin_array[atom];
            y_temp[1] = atoms::y_spin_array[atom];
            y_temp[2] = atoms::z_spin_array[atom];
            y_temp[3] = atoms::x_v_array[atom];
            y_temp[4] = atoms::y_v_array[atom];
            y_temp[5] = atoms::z_v_array[atom];
            y_temp[6] = atoms::x_w_array[atom];
            y_temp[7] = atoms::y_w_array[atom];
            y_temp[8] = atoms::z_w_array[atom];

            // Store local field in H
            std::vector<double> H = {
                atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
            };

            // Calculate k2
            auto k2 = spinDynamics(y_temp, H);
            k2_x_array[atom] = k2[0];
            k2_y_array[atom] = k2[1];
            k2_z_array[atom] = k2[2];
            k2_vx_array[atom] = k2[3];
            k2_vy_array[atom] = k2[4];
            k2_vz_array[atom] = k2[5];
            k2_wx_array[atom] = k2[6];
            k2_wy_array[atom] = k2[7];
            k2_wz_array[atom] = k2[8];

            // Update temporary state for k3
            for (size_t i = 0; i < y.size(); ++i) {
                y_temp[i] = y[i] + 0.5 * mp::dt * k2[i];
            }

            // Update arrays with temporary values
            atoms::x_spin_array[atom] = y_temp[0];
            atoms::y_spin_array[atom] = y_temp[1];
            atoms::z_spin_array[atom] = y_temp[2];
            atoms::x_v_array[atom] = y_temp[3];
            atoms::y_v_array[atom] = y_temp[4];
            atoms::z_v_array[atom] = y_temp[5];
            atoms::x_w_array[atom] = y_temp[6];
            atoms::y_w_array[atom] = y_temp[7];
            atoms::z_w_array[atom] = y_temp[8];
        }

        // Recalculate spin fields after k2 step
        calculate_spin_fields(0, num_atoms);

        // Third RK4 Step (k3)
        for (int atom = 0; atom < num_atoms; ++atom) {
            // Prepare state y_temp from previous step
            y_temp[0] = atoms::x_spin_array[atom];
            y_temp[1] = atoms::y_spin_array[atom];
            y_temp[2] = atoms::z_spin_array[atom];
            y_temp[3] = atoms::x_v_array[atom];
            y_temp[4] = atoms::y_v_array[atom];
            y_temp[5] = atoms::z_v_array[atom];
            y_temp[6] = atoms::x_w_array[atom];
            y_temp[7] = atoms::y_w_array[atom];
            y_temp[8] = atoms::z_w_array[atom];

            // Store local field in H
            std::vector<double> H = {
                atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
            };

            // Calculate k3
            auto k3 = spinDynamics(y_temp, H);
            k3_x_array[atom] = k3[0];
            k3_y_array[atom] = k3[1];
            k3_z_array[atom] = k3[2];
            k3_vx_array[atom] = k3[3];
            k3_vy_array[atom] = k3[4];
            k3_vz_array[atom] = k3[5];
            k3_wx_array[atom] = k3[6];
            k3_wy_array[atom] = k3[7];
            k3_wz_array[atom] = k3[8];

            // Update temporary state for k4
            for (size_t i = 0; i < y.size(); ++i) {
                y_temp[i] = y[i] + mp::dt * k3[i];
            }

            // Update arrays with temporary values
            atoms::x_spin_array[atom] = y_temp[0];
            atoms::y_spin_array[atom] = y_temp[1];
            atoms::z_spin_array[atom] = y_temp[2];
            atoms::x_v_array[atom] = y_temp[3];
            atoms::y_v_array[atom] = y_temp[4];
            atoms::z_v_array[atom] = y_temp[5];
            atoms::x_w_array[atom] = y_temp[6];
            atoms::y_w_array[atom] = y_temp[7];
            atoms::z_w_array[atom] = y_temp[8];
        }

        // Recalculate spin fields after k3 step
        calculate_spin_fields(0, num_atoms);

        // Fourth RK4 Step (k4)
        for (int atom = 0; atom < num_atoms; ++atom) {
            // Prepare state y_temp from previous step
            y_temp[0] = atoms::x_spin_array[atom];
            y_temp[1] = atoms::y_spin_array[atom];
            y_temp[2] = atoms::z_spin_array[atom];
            y_temp[3] = atoms::x_v_array[atom];
            y_temp[4] = atoms::y_v_array[atom];
            y_temp[5] = atoms::z_v_array[atom];
            y_temp[6] = atoms::x_w_array[atom];
            y_temp[7] = atoms::y_w_array[atom];
            y_temp[8] = atoms::z_w_array[atom];

            // Store local field in H
            std::vector<double> H = {
                atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
            };

            // Calculate k4
            auto k4 = spinDynamics(y_temp, H);
            k4_x_array[atom] = k4[0];
            k4_y_array[atom] = k4[1];
            k4_z_array[atom] = k4[2];
            k4_vx_array[atom] = k4[3];
            k4_vy_array[atom] = k4[4];
            k4_vz_array[atom] = k4[5];
            k4_wx_array[atom] = k4[6];
            k4_wy_array[atom] = k4[7];
            k4_wz_array[atom] = k4[8];
        }

        // Combine steps to get y_next
        for (int atom = 0; atom < num_atoms; ++atom) {
            std::vector<double> y_next(9);
            y_next[0] = x_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_x_array[atom] + 2 * k2_x_array[atom] + 2 * k3_x_array[atom] + k4_x_array[atom]);
            y_next[1] = y_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_y_array[atom] + 2 * k2_y_array[atom] + 2 * k3_y_array[atom] + k4_y_array[atom]);
            y_next[2] = z_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_z_array[atom] + 2 * k2_z_array[atom] + 2 * k3_z_array[atom] + k4_z_array[atom]);
            y_next[3] = v_x_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_vx_array[atom] + 2 * k2_vx_array[atom] + 2 * k3_vx_array[atom] + k4_vx_array[atom]);
            y_next[4] = v_y_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_vy_array[atom] + 2 * k2_vy_array[atom] + 2 * k3_vy_array[atom] + k4_vy_array[atom]);
            y_next[5] = v_z_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_vz_array[atom] + 2 * k2_vz_array[atom] + 2 * k3_vz_array[atom] + k4_vz_array[atom]);
            y_next[6] = w_x_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_wx_array[atom] + 2 * k2_wx_array[atom] + 2 * k3_wx_array[atom] + k4_wx_array[atom]);
            y_next[7] = w_y_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_wy_array[atom] + 2 * k2_wy_array[atom] + 2 * k3_wy_array[atom] + k4_wy_array[atom]);
            y_next[8] = w_z_initial_spin_array[atom] + (mp::dt / 6.0) * (k1_wz_array[atom] + 2 * k2_wz_array[atom] + 2 * k3_wz_array[atom] + k4_wz_array[atom]);

            // Normalize the spin length
            double S_magnitude = std::sqrt(y_next[0] * y_next[0] + y_next[1] * y_next[1] + y_next[2] * y_next[2]);
            y_next[0] /= S_magnitude;
            y_next[1] /= S_magnitude;
            y_next[2] /= S_magnitude;

            // Update atom arrays with new values
            atoms::x_spin_array[atom] = y_next[0];
            atoms::y_spin_array[atom] = y_next[1];
            atoms::z_spin_array[atom] = y_next[2];
            atoms::x_v_array[atom] = y_next[3];
            atoms::y_v_array[atom] = y_next[4];
            atoms::z_v_array[atom] = y_next[5];
            atoms::x_w_array[atom] = y_next[6];
            atoms::y_w_array[atom] = y_next[7];
            atoms::z_w_array[atom] = y_next[8];
        }

        return EXIT_SUCCESS;
    }

    // Define LLG RK4 Integrator (CUDA)
    int LLG_RK4_cuda() {
        // Add CUDA integration logic here
        return EXIT_SUCCESS;
    }
}
