//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2026. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//
//   FDT monitor (quantum:fdt-monitor): records, for a set of probe atoms, the
//   bath field b the integrator used and the spin velocity ds/dt, so that the
//   power the noise injects, P_in(w) = mu Re C_{b sdot}(w), and the power the
//   Gilbert damping removes, P_diss(w) = (alpha mu/gamma) C_{sdot sdot}(w), can
//   be compared band by band (adQTB criterion for the LLG equation). For every
//   spin dE/dt = mu b.sdot - (alpha mu/gamma)|sdot|^2 exactly; in the steady
//   state the two integrals agree, the spectra agree band by band only when the
//   system satisfies the FDT with the bath's target statistics.
//
//   Works with llg-heun on the CPU (serial) and the GPU, for the stock white
//   thermal field and for every quantum-module bath (b is read from the
//   thermal field the step used).
//
//   Output: <file>.hdr (text: sizes, units, one line per probe atom) and
//   <file> (binary, native float64): per block of n_s steps and per probe
//   atom 8 numbers
//      bx by bz      block mean of the bath field (T)
//      vx vy vz      block mean of ds/dt = (s_end - s_start)/(n_s dt) (1/s)
//      p_bs          block mean of b.sdot over the steps (exact, per step) (T/s)
//      p_ss          block mean of |sdot|^2 over the steps (1/s^2)
//   The block means act as an anti-aliasing filter: both b and sdot see the
//   same box-car, so ratios of the spectra are not affected by it.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <random>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      bool fdt_monitor = false;                        // record the FDT monitor
      std::string fdt_filename = "fdt_monitor.bin";    // binary output
      int fdt_n_probe = 256;                           // number of probe atoms
      int fdt_sample_steps = 10;                       // steps per block (n_s)
      int64_t fdt_start_step = -1;                     // first recorded step (-1: sim:equilibration-time)

      namespace fdt{

         bool initialised = false;
         bool active = false;                          // false: disabled after a setup problem
         std::vector<int> probe;                       // probe atom indices
         std::vector<double> s_old;                    // spin at the start of the step [3*np]
         std::vector<double> s_blk;                    // spin at the start of the block [3*np]
         std::vector<double> acc;                      // accumulators [8*np]: b(3), -, -, -, p_bs, p_ss
         std::vector<double> b_now;                    // this step's bath field [3*np]
         std::vector<double> s_now;                    // spin after the step [3*np]
         int n_in_block = 0;
         uint64_t n_blocks = 0;
         std::ofstream out;

         //---------------------------------------------------------------------
         // probe atoms: fixed-seed random choice, sorted (a regular stride can
         // be commensurate with the unit cell and miss sublattices)
         //---------------------------------------------------------------------
         void choose_probes(){
            const int n = atoms::num_atoms;
            const int np = std::min(fdt_n_probe, n);
            std::vector<int> all(n);
            for(int i = 0; i < n; i++) all[i] = i;
            std::mt19937 rng(20261112);
            std::shuffle(all.begin(), all.end(), rng);
            probe.assign(all.begin(), all.begin() + np);
            std::sort(probe.begin(), probe.end());
         }

         void write_header(){
            std::ofstream h((fdt_filename + ".hdr").c_str());
            const int np = probe.size();
            h << "# VAMPIRE FDT monitor header\n"
              << "n_probe " << np << "\n"
              << "sample_steps " << fdt_sample_steps << "\n"
              << "dt_SI " << mp::dt_SI << "\n"
              << "start_step " << fdt_start_step << "\n"
              << "temperature " << sim::temperature << "\n"
              << "gpu " << (gpu::acceleration ? 1 : 0) << "\n"
              << "record bx by bz vx vy vz p_bs p_ss (float64, per block, probe-major inside a block)\n"
              << "# probe: index atom material mu_s_SI(J/T) alpha gamma_SI(rad/s/T)\n";
            for(int p = 0; p < np; p++){
               const int a = probe[p];
               const int m = atoms::type_array[a];
               h << "probe " << p << " " << a << " " << m << " " << mp::material[m].mu_s_SI << " "
                 << mp::material[m].alpha << " " << mp::material[m].gamma_rel*constants::gamma_SI << "\n";
            }
         }

         void initialise(){
            initialised = true;
            if(vmpi::num_processors > 1){
               zlog << zTs() << "Warning: quantum:fdt-monitor runs in serial and GPU builds only; disabled." << std::endl;
               return;
            }
            if(fdt_start_step < 0) fdt_start_step = sim::equilibration_time;
            if(gpu::acceleration && (dipole::activated || sim::enable_fmr)){
               zlog << zTs() << "Warning: quantum:fdt-monitor on the GPU reads the bath field from the external field "
                    << "and cannot separate dipole or FMR fields; disabled." << std::endl;
               return;
            }
            choose_probes();
            const int np = probe.size();
            s_old.assign(3*np, 0.0); s_blk.assign(3*np, 0.0); acc.assign(8*np, 0.0);
            b_now.assign(3*np, 0.0); s_now.assign(3*np, 0.0);
            out.open(fdt_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            if(!out.is_open()){
               zlog << zTs() << "Warning: cannot open " << fdt_filename << "; FDT monitor disabled." << std::endl;
               return;
            }
            write_header();
            active = true;
            zlog << zTs() << "FDT monitor: " << np << " probe atoms, " << fdt_sample_steps << " steps per block, from step "
                 << fdt_start_step << ", output " << fdt_filename << std::endl;
         }

         // spins (and on the GPU the bath field) of the probe atoms
         void read_spins(std::vector<double>& s){
            if(gpu::acceleration){ gpu::fdt_gather(probe, s, b_now); return; }
            const int np = probe.size();
            for(int p = 0; p < np; p++){
               const int a = probe[p];
               s[3*p] = atoms::x_spin_array[a]; s[3*p+1] = atoms::y_spin_array[a]; s[3*p+2] = atoms::z_spin_array[a];
            }
         }

         void read_field(){
            const int np = probe.size();
            if(gpu::acceleration){   // gathered with the spins: external field = thermal + applied
               for(int p = 0; p < np; p++){
                  const int m = atoms::type_array[probe[p]];
                  const double h = mp::material[m].applied_field_strength;
                  for(int c = 0; c < 3; c++)
                     b_now[3*p+c] -= sim::H_vec[c]*sim::H_applied + h*mp::material[m].applied_field_unit_vector[c];
               }
               return;
            }
            for(int p = 0; p < np; p++){
               const int a = probe[p];
               b_now[3*p]   = atoms::x_thermal_field_array[a];
               b_now[3*p+1] = atoms::y_thermal_field_array[a];
               b_now[3*p+2] = atoms::z_thermal_field_array[a];
            }
         }

      } // end of fdt namespace

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // before the step: the spin the step starts from
   //---------------------------------------------------------------------------
   void fdt_begin_step(){
      using namespace internal;
      if(!fdt_monitor) return;
      if(!fdt::initialised) fdt::initialise();
      if(!fdt::active || static_cast<int64_t>(sim::time) < fdt_start_step) return;
      fdt::read_spins(fdt::s_old);   // GPU: b_now is overwritten here and refilled after the step
      if(fdt::n_in_block == 0) fdt::s_blk = fdt::s_old;
   }

   //---------------------------------------------------------------------------
   // after the step: bath field of the step and the new spin
   //---------------------------------------------------------------------------
   void fdt_end_step(){
      using namespace internal;
      if(!fdt_monitor || !fdt::active || static_cast<int64_t>(sim::time) < fdt_start_step) return;
      fdt::read_spins(fdt::s_now);   // GPU: also fills b_now with this step's field
      fdt::read_field();
      const int np = fdt::probe.size();
      const double idt = 1.0/mp::dt_SI;
      for(int p = 0; p < np; p++){
         double bs = 0.0, ss = 0.0;
         for(int c = 0; c < 3; c++){
            const double v = (fdt::s_now[3*p+c] - fdt::s_old[3*p+c])*idt;
            fdt::acc[8*p+c] += fdt::b_now[3*p+c];
            bs += fdt::b_now[3*p+c]*v;
            ss += v*v;
         }
         fdt::acc[8*p+6] += bs;
         fdt::acc[8*p+7] += ss;
      }
      fdt::n_in_block++;
      if(fdt::n_in_block < fdt_sample_steps) return;

      // block complete: means, write, reset
      const double inv = 1.0/fdt::n_in_block;
      const double iblk = 1.0/(fdt::n_in_block*mp::dt_SI);
      for(int p = 0; p < np; p++){
         double* r = &fdt::acc[8*p];
         for(int c = 0; c < 3; c++){
            r[c] *= inv;
            r[3+c] = (fdt::s_now[3*p+c] - fdt::s_blk[3*p+c])*iblk;
         }
         r[6] *= inv; r[7] *= inv;
      }
      fdt::out.write(reinterpret_cast<const char*>(fdt::acc.data()), fdt::acc.size()*sizeof(double));
      std::fill(fdt::acc.begin(), fdt::acc.end(), 0.0);
      fdt::n_in_block = 0;
      fdt::n_blocks++;
      if(fdt::n_blocks % 1000 == 0) fdt::out.flush();
   }

   void fdt_finalize(){
      if(internal::fdt::out.is_open()){
         internal::fdt::out.close();
         zlog << zTs() << "FDT monitor: " << internal::fdt::n_blocks << " blocks written" << std::endl;
      }
   }

} // end of quantum namespace
