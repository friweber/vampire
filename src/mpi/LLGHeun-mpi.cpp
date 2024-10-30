#ifdef MPICF
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include <iomanip>

#include <cmath>
#include <vector>

int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);
int set_LLG();

namespace sim {

int LLG_Heun_mpi() {
    // Check for error
    if(err::check==true){std::cout << "LLG_Heun_mpi has been called" << std::endl;}

	using namespace LLG_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGinit();

    const int num_atoms = atoms::num_atoms;
    const int pre_comm_si = 0;
    const int pre_comm_ei = vmpi::num_core_atoms;
    const int post_comm_si = vmpi::num_core_atoms;
    const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

    // Constants for Non-Markovian LLG
    const double A = 152.e5;
    const double omega0 = 282.5;
    const double Gamma = 166.7;
    const double gyromagnetic_ratio = 1.;
    calculate_external_fields(0, num_atoms);
    // Initiate halo swap
    vmpi::mpi_init_halo_swap();

    // Store initial spin, V, and W positions
    for(int atom=pre_comm_si; atom<post_comm_ei; atom++){
        x_initial_spin_array[atom] = atoms::x_spin_array[atom];
        y_initial_spin_array[atom] = atoms::y_spin_array[atom];
        z_initial_spin_array[atom] = atoms::z_spin_array[atom];
        v_x_initial_spin_array[atom] = atoms::x_v_array[atom];
        v_y_initial_spin_array[atom] = atoms::y_v_array[atom];
        v_z_initial_spin_array[atom] = atoms::z_v_array[atom];
        w_x_initial_spin_array[atom] = atoms::x_w_array[atom];
        w_y_initial_spin_array[atom] = atoms::y_w_array[atom];
        w_z_initial_spin_array[atom] = atoms::z_w_array[atom];
    }

    // Calculate fields (core)
    calculate_spin_fields(pre_comm_si, pre_comm_ei);

    // Calculate Euler Step (Core)
    for(int atom=pre_comm_si; atom<pre_comm_ei; atom++){
        const int imaterial = atoms::type_array[atom];

        const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
        const double V[3] = {atoms::x_v_array[atom], atoms::y_v_array[atom], atoms::z_v_array[atom]};
        const double W[3] = {atoms::x_w_array[atom], atoms::y_w_array[atom], atoms::z_w_array[atom]};
        const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                             atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                             atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

        // Calculate dS/dt
        double dS[3];
        dS[0] = gyromagnetic_ratio * (S[1]*(H[2]+V[2]) - S[2]*(H[1]+V[1]));
        dS[1] = gyromagnetic_ratio * (S[2]*(H[0]+V[0]) - S[0]*(H[2]+V[2]));
        dS[2] = gyromagnetic_ratio * (S[0]*(H[1]+V[1]) - S[1]*(H[0]+V[0]));

        // Calculate dV/dt and dW/dt
        double dV[3], dW[3];
        for(int i=0; i<3; i++){
            dV[i] = W[i];
            dW[i] = -omega0*omega0 * V[i] - Gamma * W[i] + A * gyromagnetic_ratio * S[i];
        }

        // Store dS, dV, dW in euler arrays
        x_euler_array[atom] = dS[0];
        y_euler_array[atom] = dS[1];
        z_euler_array[atom] = dS[2];

        w_x_euler_array[atom] = dW[0];
        w_y_euler_array[atom] = dW[1];
        w_z_euler_array[atom] = dW[2];

        v_x_euler_array[atom] = dV[0];
        v_y_euler_array[atom] = dV[1];
        v_z_euler_array[atom] = dV[2];

        // Calculate Euler Step
        double S_new[3], V_new[3], W_new[3];
        for(int i=0; i<3; i++){
            S_new[i] = S[i] + dS[i] * material_parameters::dt;
            V_new[i] = V[i] + dV[i] * material_parameters::dt;
            W_new[i] = W[i] + dW[i] * material_parameters::dt;
        }

        // Normalize Spin Length
        double mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
        for(int i=0; i<3; i++) S_new[i] *= mod_S;

        // Store new values
        x_spin_storage_array[atom] = S_new[0];
        y_spin_storage_array[atom] = S_new[1];
        z_spin_storage_array[atom]= S_new[2];
        v_x_spin_storage_array[atom] = V_new[0];
        v_y_spin_storage_array[atom] = V_new[1];
        v_z_spin_storage_array[atom] = V_new[2];
        w_x_spin_storage_array[atom] = W_new[0];
        w_y_spin_storage_array[atom] = W_new[1];
        w_z_spin_storage_array[atom] = W_new[2];
    }

    // Complete halo swap
    vmpi::mpi_complete_halo_swap();

    // Calculate fields (boundary)
    calculate_spin_fields(post_comm_si, post_comm_ei);

    // Calculate Euler Step (boundary)
    for(int atom=post_comm_si; atom<post_comm_ei; atom++){
        const int imaterial = atoms::type_array[atom];

        const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
        const double V[3] = {atoms::x_v_array[atom], atoms::y_v_array[atom], atoms::z_v_array[atom]};
        const double W[3] = {atoms::x_w_array[atom], atoms::y_w_array[atom], atoms::z_w_array[atom]};
        const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                             atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                             atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

        // Calculate dS/dt
        double dS[3];
        dS[0] = gyromagnetic_ratio * (S[1]*(H[2]+V[2]) - S[2]*(H[1]+V[1]));
        dS[1] = gyromagnetic_ratio * (S[2]*(H[0]+V[0]) - S[0]*(H[2]+V[2]));
        dS[2] = gyromagnetic_ratio * (S[0]*(H[1]+V[1]) - S[1]*(H[0]+V[0]));

        // Calculate dV/dt and dW/dt
        double dV[3], dW[3];
        for(int i=0; i<3; i++){
            dV[i] = W[i];
            dW[i] = -omega0*omega0 * V[i] - Gamma * W[i] + A * gyromagnetic_ratio * S[i];
        }

        // Store dS, dV, dW in euler arrays
        x_euler_array[atom] = dS[0];
        y_euler_array[atom] = dS[1];
        z_euler_array[atom] = dS[2];

        w_x_euler_array[atom] = dW[0];
        w_y_euler_array[atom] = dW[1];
        w_z_euler_array[atom] = dW[2];

        v_x_euler_array[atom] = dV[0];
        v_y_euler_array[atom] = dV[1];
        v_z_euler_array[atom] = dV[2];

        // Calculate Euler Step
        double S_new[3], V_new[3], W_new[3];
        for(int i=0; i<3; i++){
            S_new[i] = S[i] + dS[i] * material_parameters::dt;
            V_new[i] = V[i] + dV[i] * material_parameters::dt;
            W_new[i] = W[i] + dW[i] * material_parameters::dt;
        }

        // Normalize Spin Length
        double mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
        for(int i=0; i<3; i++) S_new[i] *= mod_S;

        // Store new values
        x_spin_storage_array[atom] = S_new[0];
        y_spin_storage_array[atom] = S_new[1];
        z_spin_storage_array[atom]= S_new[2];
        v_x_spin_storage_array[atom] = V_new[0];
        v_y_spin_storage_array[atom] = V_new[1];
        v_z_spin_storage_array[atom] = V_new[2];
        w_x_spin_storage_array[atom] = W_new[0];
        w_y_spin_storage_array[atom] = W_new[1];
        w_z_spin_storage_array[atom] = W_new[2];
    }

    // Copy new spins to spin array (all)
    for(int atom=pre_comm_si; atom<post_comm_ei; atom++){
        atoms::x_spin_array[atom] = x_spin_storage_array[atom];
        atoms::y_spin_array[atom] = y_spin_storage_array[atom];
        atoms::z_spin_array[atom] = z_spin_storage_array[atom];
        atoms::x_v_array[atom] = v_x_spin_storage_array[atom];
        atoms::y_v_array[atom] = v_y_spin_storage_array[atom];
        atoms::z_v_array[atom] = v_z_spin_storage_array[atom];
        atoms::x_w_array[atom] = w_x_spin_storage_array[atom];
        atoms::y_w_array[atom] = w_y_spin_storage_array[atom];
        atoms::z_w_array[atom] =  w_z_spin_storage_array[atom];
    }

    // Initiate second halo swap
    vmpi::mpi_init_halo_swap();

    // Recalculate spin dependent fields (core)
    calculate_spin_fields(pre_comm_si, pre_comm_ei);

    // Calculate Heun Gradients (core)
    for(int atom=pre_comm_si; atom<pre_comm_ei; atom++){
        const int imaterial = atoms::type_array[atom];


        const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
        const double V[3] = {atoms::x_v_array[atom], atoms::y_v_array[atom], atoms::z_v_array[atom]};
        const double W[3] = {atoms::x_w_array[atom], atoms::y_w_array[atom], atoms::z_w_array[atom]};
        const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                             atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                             atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

        // Calculate dS/dt, dV/dt, dW/dt (same as Euler step)
        double dS[3], dV[3], dW[3];
        dS[0] = gyromagnetic_ratio * (S[1]*(H[2]+V[2]) - S[2]*(H[1]+V[1]));
        dS[1] = gyromagnetic_ratio * (S[2]*(H[0]+V[0]) - S[0]*(H[2]+V[2]));
        dS[2] = gyromagnetic_ratio * (S[0]*(H[1]+V[1]) - S[1]*(H[0]+V[0]));
        for(int i=0; i<3; i++){
            dV[i] = W[i];
            dW[i] = -omega0*omega0 * V[i] - Gamma * W[i] + A * gyromagnetic_ratio * S[i];
        }

        // Store dS in heun array
        x_heun_array[atom] = dS[0];
        y_heun_array[atom] = dS[1];
        z_heun_array[atom] = dS[2];

        v_x_heun_array[atom] = dV[0];
        v_y_heun_array[atom] = dV[1];
        v_z_heun_array[atom] = dV[2];

        w_x_heun_array[atom] = dW[0];
        w_y_heun_array[atom] = dW[1];
        w_z_heun_array[atom] = dW[2];
    }


    // Complete second halo swap
    vmpi::mpi_complete_halo_swap();

    //----------------------------------------
	// Calculate Heun Gradients (boundary)
	//----------------------------------------


    
	for(int atom=post_comm_si;atom<post_comm_ei;atom++){
        const int imaterial = atoms::type_array[atom];


        const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
        const double V[3] = {atoms::x_v_array[atom], atoms::y_v_array[atom], atoms::z_v_array[atom]};
        const double W[3] = {atoms::x_w_array[atom], atoms::y_w_array[atom], atoms::z_w_array[atom]};
        const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                             atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                             atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};

        // Calculate dS/dt, dV/dt, dW/dt (same as Euler step)
        double dS[3], dV[3], dW[3];
        dS[0] = gyromagnetic_ratio * (S[1]*(H[2]+V[2]) - S[2]*(H[1]+V[1]));
        dS[1] = gyromagnetic_ratio * (S[2]*(H[0]+V[0]) - S[0]*(H[2]+V[2]));
        dS[2] = gyromagnetic_ratio * (S[0]*(H[1]+V[1]) - S[1]*(H[0]+V[0]));
        for(int i=0; i<3; i++){
            dV[i] = W[i];
            dW[i] = -omega0*omega0 * V[i] - Gamma * W[i] + A * gyromagnetic_ratio * S[i];
        }

        // Store dS in heun array
        x_heun_array[atom] = dS[0];
        y_heun_array[atom] = dS[1];
        z_heun_array[atom] = dS[2];

        v_x_heun_array[atom] = dV[0];
        v_y_heun_array[atom] = dV[1];
        v_z_heun_array[atom] = dV[2];

        w_x_heun_array[atom] = dW[0];
        w_y_heun_array[atom] = dW[1];
        w_z_heun_array[atom] = dW[2];
    }


    // Calculate Heun Step
    for(int atom=pre_comm_si; atom<post_comm_ei; atom++){
        double S_new[3], V_new[3], W_new[3];
        
        // Calculate new S
        S_new[0] = x_initial_spin_array[atom] + material_parameters::half_dt * (x_euler_array[atom] + x_heun_array[atom]);
        S_new[1] = y_initial_spin_array[atom] + material_parameters::half_dt * (y_euler_array[atom] + y_heun_array[atom]);
        S_new[2] = z_initial_spin_array[atom] + material_parameters::half_dt * (z_euler_array[atom] + z_heun_array[atom]);

        
        V_new[0] = v_x_initial_spin_array[atom] + material_parameters::half_dt * (v_x_euler_array[atom] + v_x_heun_array[atom]);
        V_new[1] = v_y_initial_spin_array[atom] + material_parameters::half_dt * (v_y_euler_array[atom] + v_y_heun_array[atom]);
        V_new[2] = v_z_initial_spin_array[atom] + material_parameters::half_dt * (v_z_euler_array[atom] + v_z_heun_array[atom]);
       
       
        W_new[0] = w_x_initial_spin_array[atom] + material_parameters::half_dt * (w_x_euler_array[atom] + w_x_heun_array[atom]);
        W_new[1] = w_y_initial_spin_array[atom] + material_parameters::half_dt * (w_y_euler_array[atom] + w_y_heun_array[atom]);
        W_new[2] = w_z_initial_spin_array[atom] + material_parameters::half_dt * (w_z_euler_array[atom] + w_z_heun_array[atom]);

        // Normalize Spin Length
        double mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
        S_new[0] *= mod_S;
        S_new[1] *= mod_S;
        S_new[2] *= mod_S;

        // Copy new spins to spin array
        atoms::x_spin_array[atom] = S_new[0];
        atoms::y_spin_array[atom] = S_new[1];
        atoms::z_spin_array[atom] = S_new[2];
        atoms::x_v_array[atom] = V_new[0];
        atoms::y_v_array[atom] = V_new[1];
        atoms::z_v_array[atom] = V_new[2];
        atoms::x_w_array[atom] = W_new[0];
        atoms::y_w_array[atom] = W_new[1];
        atoms::z_w_array[atom] = W_new[2];
    }

    // Swap timers compute -> wait
    vmpi::TotalComputeTime += vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

    // Wait for other processors
    vmpi::barrier();

    // Swap timers wait -> compute
    vmpi::TotalWaitTime += vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);


    // const int atom = 0;
    // std::cout << std::setprecision(17);
    // std::cout << "Initial spin positions for atom 0:" << std::endl;
    // std::cout << "Sx: " << x_initial_spin_array[atom] << ", Sy: " << y_initial_spin_array[atom] << ", Sz: " << z_initial_spin_array[atom] << std::endl;
    // std::cout << "Vx: " << v_x_initial_spin_array[atom] << ", Vy: " << v_y_initial_spin_array[atom] << ", Vz: " << v_z_initial_spin_array[atom] << std::endl;
    // std::cout << "Wx: " << w_x_initial_spin_array[atom] << ", Wy: " << w_y_initial_spin_array[atom] << ", Wz: " << w_z_initial_spin_array[atom] << std::endl;

    // std::cout << "Field values for atom 0:" << std::endl;
    // std::cout << "Hx: " << atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] << ", Hy: " << atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] << ", Hz: " << atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] << std::endl;

    // std::cout << "Euler step results for atom 0:" << std::endl;
    // std::cout << "dSx: " << x_euler_array[atom] << ", dSy: " << y_euler_array[atom] << ", dSz: " << z_euler_array[atom] << std::endl;
    // std::cout << "dVx: " << v_x_euler_array[atom] << ", dVy: " << v_y_euler_array[atom] << ", dVz: " << v_z_euler_array[atom] << std::endl;
    // std::cout << "dWx: " << w_x_euler_array[atom] << ", dWy: " << w_y_euler_array[atom] << ", dWz: " << w_z_euler_array[atom] << std::endl;

    // std::cout << "Heun step results for atom 0:" << std::endl;
    // std::cout << "dSx: " << x_heun_array[atom] << ", dSy: " << y_heun_array[atom] << ", dSz: " << z_heun_array[atom] << std::endl;
    // std::cout << "dVx: " << v_x_heun_array[atom] << ", dVy: " << v_y_heun_array[atom] << ", dVz: " << v_z_heun_array[atom] << std::endl;
    // std::cout << "dWx: " << w_x_heun_array[atom] << ", dWy: " << w_y_heun_array[atom] << ", dWz: " << w_z_heun_array[atom] << std::endl;

    // std::cout << "Final spin positions for atom 0:" << std::endl;
    // std::cout << "Sx: " << atoms::x_spin_array[atom] << ", Sy: " << atoms::y_spin_array[atom] << ", Sz: " << atoms::z_spin_array[atom] << std::endl;
    // std::cout << "Vx: " << atoms::x_v_array[atom] << ", Vy: " << atoms::y_v_array[atom] << ", Vz: " << atoms::z_v_array[atom] << std::endl;
    // std::cout << "Wx: " << atoms::x_w_array[atom] << ", Wy: " << atoms::y_w_array[atom] << ", Wz: " << atoms::z_w_array[atom] << std::endl;

    return EXIT_SUCCESS;
}

} // end of namespace sim
#endif