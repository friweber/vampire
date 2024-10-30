/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
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
int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);


namespace LLG_arrays {
    // Local arrays for LLG integration
    std::vector<double> x_euler_array;
    std::vector<double> y_euler_array;
    std::vector<double> z_euler_array;

    std::vector<double> x_heun_array;
    std::vector<double> y_heun_array;
    std::vector<double> z_heun_array;

    std::vector<double> x_spin_storage_array;
    std::vector<double> y_spin_storage_array;
    std::vector<double> z_spin_storage_array;

    std::vector<double> x_initial_spin_array;
    std::vector<double> y_initial_spin_array;
    std::vector<double> z_initial_spin_array;

    // Non-Markovian variables
    std::vector<double> w_x_euler_array;
    std::vector<double> w_y_euler_array;
    std::vector<double> w_z_euler_array;

    std::vector<double> w_x_heun_array;
    std::vector<double> w_y_heun_array;
    std::vector<double> w_z_heun_array;

    std::vector<double> w_x_spin_storage_array;
    std::vector<double> w_y_spin_storage_array;
    std::vector<double> w_z_spin_storage_array;

    std::vector<double> w_x_initial_spin_array;
    std::vector<double> w_y_initial_spin_array;
    std::vector<double> w_z_initial_spin_array;

    // Local arrays for LLG integration for components of v
    std::vector<double> v_x_euler_array;
    std::vector<double> v_y_euler_array;
    std::vector<double> v_z_euler_array;

    std::vector<double> v_x_heun_array;
    std::vector<double> v_y_heun_array;
    std::vector<double> v_z_heun_array;

    std::vector<double> v_x_spin_storage_array;
    std::vector<double> v_y_spin_storage_array;
    std::vector<double> v_z_spin_storage_array;

    std::vector<double> v_x_initial_spin_array;
    std::vector<double> v_y_initial_spin_array;
    std::vector<double> v_z_initial_spin_array;

    bool LLG_set = false; ///< Flag to define state of LLG arrays (initialized/uninitialized)
}

namespace sim {
    int LLGinit() {
        // Check calling of routine if error checking is activated
        if (err::check == true) { std::cout << "sim::LLG_init has been called" << std::endl; }

        using namespace LLG_arrays;

        std::cout << "Non Markovian LLG has been called" << std::endl;


        x_spin_storage_array.resize(atoms::num_atoms, 0.0);
        y_spin_storage_array.resize(atoms::num_atoms, 0.0);
        z_spin_storage_array.resize(atoms::num_atoms, 0.0);

        x_initial_spin_array.resize(atoms::num_atoms, 0.0);
        y_initial_spin_array.resize(atoms::num_atoms, 0.0);
        z_initial_spin_array.resize(atoms::num_atoms, 0.0);

        x_euler_array.resize(atoms::num_atoms, 0.0);
        y_euler_array.resize(atoms::num_atoms, 0.0);
        z_euler_array.resize(atoms::num_atoms, 0.0);

        x_heun_array.resize(atoms::num_atoms, 0.0);
        y_heun_array.resize(atoms::num_atoms, 0.0);
        z_heun_array.resize(atoms::num_atoms, 0.0);

        // Storage arrays for components of w
        w_x_spin_storage_array.resize(atoms::num_atoms, 0.0);
        w_y_spin_storage_array.resize(atoms::num_atoms, 0.0);
        w_z_spin_storage_array.resize(atoms::num_atoms, 0.0);

        // Storage arrays for components of v
        v_x_spin_storage_array.resize(atoms::num_atoms, 0.0);
        v_y_spin_storage_array.resize(atoms::num_atoms, 0.0);
        v_z_spin_storage_array.resize(atoms::num_atoms, 0.0);

        // Initial spin arrays for components of w
        w_x_initial_spin_array.resize(atoms::num_atoms, 0.0);
        w_y_initial_spin_array.resize(atoms::num_atoms, 0.0);
        w_z_initial_spin_array.resize(atoms::num_atoms, 0.0);

        // Initial spin arrays for components of v
        v_x_initial_spin_array.resize(atoms::num_atoms, 0.0);
        v_y_initial_spin_array.resize(atoms::num_atoms, 0.0);
        v_z_initial_spin_array.resize(atoms::num_atoms, 0.0);

        // Euler method arrays for components of w
        w_x_euler_array.resize(atoms::num_atoms, 0.0);
        w_y_euler_array.resize(atoms::num_atoms, 0.0);
        w_z_euler_array.resize(atoms::num_atoms, 0.0);

        // Euler method arrays for components of v
        v_x_euler_array.resize(atoms::num_atoms, 0.0);
        v_y_euler_array.resize(atoms::num_atoms, 0.0);
        v_z_euler_array.resize(atoms::num_atoms, 0.0);

        // Heun method arrays for components of w
        w_x_heun_array.resize(atoms::num_atoms, 0.0);
        w_y_heun_array.resize(atoms::num_atoms, 0.0);
        w_z_heun_array.resize(atoms::num_atoms, 0.0);

        // Heun method arrays for components of v
        v_x_heun_array.resize(atoms::num_atoms, 0.0);
        v_y_heun_array.resize(atoms::num_atoms, 0.0);
        v_z_heun_array.resize(atoms::num_atoms, 0.0);

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
    std::vector<double> spinDynamics(const std::vector<double>& y,const std::vector<double>& b) {
        const  double A = 152.e5;
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

    std::vector<double> Euler_step(
    std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)> ode,
    const std::vector<double>& y, 
    const std::vector<double>& b_start, 
    const double& dt,
    std::vector<double>& k1_out) {

    k1_out = ode(y, b_start);
    std::vector<double> y_pred(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y_pred[i] = y[i] + dt * k1_out[i];
    }

    return y_pred;
    }

    std::vector<double> Heun_step(
    std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)> ode,
    const std::vector<double>& y, 
    const std::vector<double>& y_pred,
    const std::vector<double>& b_end, 
    const double& dt,
    const std::vector<double>& k1) {

    std::vector<double> k2 = ode(y_pred, b_end);
    std::vector<double> y_next(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y_next[i] = y[i] + 0.5 * dt * (k1[i] + k2[i]);
    }

    return y_next;
    }


    int LLG_Heun() {
    using namespace LLG_arrays;

    if (LLG_set == false) sim::LLGinit();

    const int num_atoms = atoms::num_atoms;
    std::vector<double> y(9);



    calculate_spin_fields(0, num_atoms);
    calculate_external_fields(0, num_atoms);

    int atom = 0; // First atom

    const double S_x = atoms::x_spin_array[atom];
    const double S_y = atoms::y_spin_array[atom];
    const double S_z = atoms::z_spin_array[atom];

    const double V_x = atoms::x_v_array[atom];
    const double V_y = atoms::y_v_array[atom];
    const double V_z = atoms::z_v_array[atom];

    const double W_x = atoms::x_w_array[atom];
    const double W_y = atoms::y_w_array[atom];
    const double W_z = atoms::z_w_array[atom];

    const double H_ext_x = atoms::x_total_external_field_array[atom];
    const double H_ext_y = atoms::y_total_external_field_array[atom];
    const double H_ext_z = atoms::z_total_external_field_array[atom];

    const double H_spin_x = atoms::x_total_spin_field_array[atom];
    const double H_spin_y = atoms::y_total_spin_field_array[atom];
    const double H_spin_z = atoms::z_total_spin_field_array[atom];

    std::cout << std::setprecision(17);
    std::cout << S_x << " " << S_y << " " << S_z << " ";
    std::cout << V_x << " " << V_y << " " << V_z << " ";
    std::cout << W_x << " " << W_y << " " << W_z << " ";
    std::cout << H_ext_x << " " << H_ext_y << " " << H_ext_z << " ";
    std::cout << H_spin_x << " " << H_spin_y << " " << H_spin_z << std::endl;

    int actual_realizations = atoms::noise_field.size();

    std::vector<std::vector<double>> k1_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> k2_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> k3_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> k4_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> y_pred_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> y_in_storage(num_atoms, std::vector<double>(9));

    for (int atom = 0; atom < num_atoms; atom++) {
        y_in_storage[atom][0] = atoms::x_spin_array[atom];
        y_in_storage[atom][1] = atoms::y_spin_array[atom];
        y_in_storage[atom][2] = atoms::z_spin_array[atom];
        y_in_storage[atom][3] = atoms::x_v_array[atom];
        y_in_storage[atom][4] = atoms::y_v_array[atom];
        y_in_storage[atom][5] = atoms::z_v_array[atom];
        y_in_storage[atom][6] = atoms::x_w_array[atom];
        y_in_storage[atom][7] = atoms::y_w_array[atom];
        y_in_storage[atom][8] = atoms::z_w_array[atom];
    }


    // K1 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        y[0] = atoms::x_spin_array[atom];
        y[1] = atoms::y_spin_array[atom];
        y[2] = atoms::z_spin_array[atom];
        y[3] = atoms::x_v_array[atom];
        y[4] = atoms::y_v_array[atom];
        y[5] = atoms::z_v_array[atom];
        y[6] = atoms::x_w_array[atom];
        y[7] = atoms::y_w_array[atom];
        y[8] = atoms::z_w_array[atom];

        std::vector<double> H = {
            atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
            atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
            atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
        };


        k1_storage[atom] = spinDynamics(y, H);

        for (size_t i = 0; i < y.size(); ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + 0.5 * mp::dt * k1_storage[atom][i];
        }


        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];


    }

    calculate_spin_fields(0, num_atoms);

    // K2 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        y[0] = atoms::x_spin_array[atom];
        y[1] = atoms::y_spin_array[atom];
        y[2] = atoms::z_spin_array[atom];
        y[3] = atoms::x_v_array[atom];
        y[4] = atoms::y_v_array[atom];
        y[5] = atoms::z_v_array[atom];
        y[6] = atoms::x_w_array[atom];
        y[7] = atoms::y_w_array[atom];
        y[8] = atoms::z_w_array[atom];

        std::vector<double> H = {
            atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
            atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
            atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
        };


        k2_storage[atom] = spinDynamics(y, H);

        for (size_t i = 0; i < y.size(); ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + 0.5 * mp::dt * k2_storage[atom][i];
        }



        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];


    }

    calculate_spin_fields(0, num_atoms);

    // K3 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        y[0] = atoms::x_spin_array[atom];
        y[1] = atoms::y_spin_array[atom];
        y[2] = atoms::z_spin_array[atom];
        y[3] = atoms::x_v_array[atom];
        y[4] = atoms::y_v_array[atom];
        y[5] = atoms::z_v_array[atom];
        y[6] = atoms::x_w_array[atom];
        y[7] = atoms::y_w_array[atom];
        y[8] = atoms::z_w_array[atom];

        std::vector<double> H = {
            atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
            atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
            atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
        };


        k3_storage[atom] = spinDynamics(y, H);

        for (size_t i = 0; i < y.size(); ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + mp::dt * k3_storage[atom][i];
        }


        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];


    }

    calculate_spin_fields(0, num_atoms);
    
    //K4 step
    for (int atom = 0; atom < num_atoms; ++atom) {
        y[0] = atoms::x_spin_array[atom];
        y[1] = atoms::y_spin_array[atom];
        y[2] = atoms::z_spin_array[atom];
        y[3] = atoms::x_v_array[atom];
        y[4] = atoms::y_v_array[atom];
        y[5] = atoms::z_v_array[atom];
        y[6] = atoms::x_w_array[atom];
        y[7] = atoms::y_w_array[atom];
        y[8] = atoms::z_w_array[atom];

        std::vector<double> H = {
            atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
            atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
            atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]
        };


        k4_storage[atom] = spinDynamics(y, H);

        for (size_t i = 0; i < y.size(); ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + (mp::dt / 6.0) * (k1_storage[atom][i] + 2.0 * k2_storage[atom][i] + 2.0 * k3_storage[atom][i] + k4_storage[atom][i]);
        }


        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];


    }

    return EXIT_SUCCESS;
}



}
