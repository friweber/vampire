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
//   Noise generation and retrieval functions for the quantum thermostat.
//
//   This file contains:
//     - PSD():                  Power spectral density evaluation
//     - assign_unique_indices(): Map atoms to flat noise array offsets
//     - calculate_noise():      Non-windowed full-run noise generation
//     - get_noise():            Linear interpolation on coarse noise grid
//     - get_field():            Public API to retrieve noise for an atom
//     - increment_time():       Advance the fine step counter
//     - export_noise_data():    Append atom-0 noise from coarse_noise_field to file
//     - init_noise_structures(): Setup persistent FFT resources (windowed)
//     - generate_noise_window(): Generate one window of noise
//     - update_noise_if_needed(): Check and regenerate window as needed
//     - cleanup_noise_structures(): Free persistent FFT resources
//
//   Noise generation pipeline (FFT method):
//     1. Generate white noise with sigma = 1/sqrt(dt_fine)
//     2. Forward FFT to frequency domain
//     3. Multiply by sqrt(PSD) to impose Lorentzian spectral shape
//     4. Inverse FFT back to time domain
//     5. Scale by: (1/n_coarse) * (1/sqrt(S0)) * sqrt(dt_fine/dt_coarse)
//     6. Store in coarse_noise_field for later interpolation
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#ifdef FFT
#include <fftw3.h>
#endif

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// Module headers
#include "internal.hpp"

namespace quantum{
   namespace internal{

      //=====================================================================
      // Shared FFT shaping body for both calculate_noise() (full-run) and
      // generate_noise_window() (overlap-save windowed).
      //
      // Forward FFT on in_buf → multiply each frequency bin by sqrt_psd[i]
      // → inverse FFT into out_time. The caller is responsible for filling
      // in_buf (white noise, possibly with an overlap tail) and for any
      // post-processing of out_time (e.g., output scaling, extracting the
      // valid middle).
      //
      // nc = full FFT length. sqrt_psd has length nc/2+1.
      //=====================================================================
      void run_fft_pipeline(double* in_buf,
                            fftw_plan plan_fwd, fftw_complex* out_freq,
                            fftw_plan plan_bwd, double* out_time,
                            const std::vector<double>& sqrt_psd, int nc) {
         #ifdef FFT
         (void)in_buf;    // plans bind to the buffers; args kept for clarity
         (void)out_time;
         fftw_execute(plan_fwd);
         for (int i = 0; i <= nc / 2; ++i) {
            out_freq[i][0] *= sqrt_psd[i];
            out_freq[i][1] *= sqrt_psd[i];
         }
         fftw_execute(plan_bwd);
         #else
         (void)in_buf; (void)plan_fwd; (void)out_freq;
         (void)plan_bwd; (void)out_time; (void)sqrt_psd; (void)nc;
         #endif
      }

      //=====================================================================
      // Power Spectral Density
      //
      // Evaluates the Lorentzian PSD at angular frequency omega for
      // temperature T and the given material. Three modes are supported:
      //   classical:     P = 2T * A * Gamma / D
      //   quantum:       P = coth(omega/2T) * A * Gamma * omega / D
      //   quantum-no-zero: P = [coth(omega/2T) - 1] * A * Gamma * omega / D
      // where D = (omega0^2 - omega^2)^2 + Gamma^2 * omega^2
      //=====================================================================
      double PSD(const double omega, const double T, const int material) {

         const double A      = material_A_array[material];
         const double Gamma  = material_gamma_array[material];
         const double omega0 = material_omega0_array[material];

         // Lorentzian denominator: (omega0^2 - omega^2)^2 + Gamma^2 * omega^2
         double lorentzian_denom = (omega0*omega0 - omega*omega) * (omega0*omega0 - omega*omega)
                                 + Gamma*Gamma * omega*omega;
         if (lorentzian_denom < 1e-12) lorentzian_denom = 1e-12;

         // coth(x) with stabilization near x=0
         double x = (T > 1e-12) ? omega / (2.0 * T) : omega;
         double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);

         switch (internal::noise_type) {

            case internal::classical:
               return 2.0 * T * A * Gamma / lorentzian_denom;

            case internal::quantum_zero:
               if (omega > 0) return coth * A * Gamma * omega / lorentzian_denom;
               else return 2.0 * T * A * Gamma / (omega0*omega0 * omega0*omega0);

            case internal::quantum_no_zero:
               if (omega > 0) return (coth - 1.0) * A * Gamma * omega / lorentzian_denom;
               else return 2.0 * T * A * Gamma / (omega0*omega0 * omega0*omega0);

            default:
               zlog << zTs() << "Programmer error: unknown quantum noise type " << internal::noise_type << std::endl;
               std::cerr << "Programmer error: unknown quantum noise type " << internal::noise_type << std::endl;
               err::vexit();
               return 0.0;
         }
      }

      //=====================================================================
      // Assign unique indices for each atom and spatial component
      //
      // The coarse_noise_field is a flat 1D array laid out as:
      //   [atom0_x(0..nc-1), atom0_y(0..nc-1), atom0_z(0..nc-1),
      //    atom1_x(0..nc-1), atom1_y(0..nc-1), atom1_z(0..nc-1), ...]
      //
      // This function precomputes the starting offset for each atom's
      // x, y, z noise data so that get_noise() can quickly look up values.
      //=====================================================================
      void assign_unique_indices(int n_coarse, int num_atoms_local) {

         std::cout << "Assigning indices for " << num_atoms_local
                   << " local atoms with " << n_coarse << " coarse steps." << std::endl;

         atom_idx_x.resize(num_atoms_local);
         atom_idx_y.resize(num_atoms_local);
         atom_idx_z.resize(num_atoms_local);

         for (int atom = 0; atom < num_atoms_local; atom++) {
            const size_t a = static_cast<size_t>(atom);
            const size_t nc = static_cast<size_t>(n_coarse);
            atom_idx_x[atom] = 3 * a * nc;
            atom_idx_y[atom] = 3 * a * nc + nc;
            atom_idx_z[atom] = 3 * a * nc + 2 * nc;
         }
      }

      //=====================================================================
      // Non-windowed noise generation (full simulation at once)
      //
      // Generates colored noise for all realizations on a coarse time grid,
      // using a temporary FFT. The noise_field vector is resized and filled.
      //
      // Noise scaling chain:
      //   white noise sigma = 1/sqrt(dt_fine)
      //   -> FFT -> multiply by sqrt(PSD) -> iFFT
      //   -> scale by (1/n_coarse) * (1/sqrt(S0)) * sqrt(dt_fine/dt_coarse)
      //=====================================================================
      void calculate_noise(int realizations, double dt_fine,
                           int M, double T, int n_coarse_total,
                           std::vector<double>& noise_field) {
         #ifdef FFT

         const double dt_coarse = dt_fine * M;
         const double S0 = material_S0_array[0];
         const double inv_sqrt_S0 = (S0 > 0.0) ? 1.0 / std::sqrt(S0) : 1.0;

         // Temporary FFTW arrays for the full coarse grid
         double* __restrict in       = (double*)fftw_malloc(sizeof(double) * n_coarse_total);
         fftw_complex* __restrict out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n_coarse_total/2 + 1));
         double* __restrict result    = (double*)fftw_malloc(sizeof(double) * n_coarse_total);

         fftw_plan forward  = fftw_plan_dft_r2c_1d(n_coarse_total, in, out, FFTW_MEASURE);
         fftw_plan backward = fftw_plan_dft_c2r_1d(n_coarse_total, out, result, FFTW_MEASURE);

         const double norm_factor = 1.0 / n_coarse_total;
         const double sigma = 1.0 / std::sqrt(dt_fine);
         const double scale = std::sqrt(dt_fine / dt_coarse);

         // Allocate output vector
         try {
            const size_t total_elements = static_cast<size_t>(realizations) * static_cast<size_t>(n_coarse_total);
            noise_field.resize(total_elements);
            std::cout << "  Noise field allocated: " << total_elements
                      << " elements (" << (total_elements*sizeof(double))/(1024*1024) << " MB)" << std::endl;
         } catch (const std::exception& e) {
            std::cerr << "Error allocating noise field: " << e.what() << std::endl;
            err::vexit();
         }

         // Precompute sqrt(PSD) on the coarse frequency grid
         std::vector<double> sqrt_PSD_coarse(n_coarse_total/2 + 1);
         const double df_coarse = 1.0 / (n_coarse_total * dt_coarse);
         for (int i = 0; i <= n_coarse_total/2; ++i) {
            double omega = 2.0 * M_PI * i * df_coarse;
            sqrt_PSD_coarse[i] = std::sqrt(PSD(omega, T, 0));
         }

         std::cout << "  Generating noise for " << realizations << " realizations"
                   << " (" << n_coarse_total << " coarse steps)..." << std::endl;

         // Progress tracking
         const int bar_width = 50;
         int last_printed_percent = -1;

         for (int r = 0; r < realizations; ++r) {

            // 1. Generate white noise with sigma = 1/sqrt(dt_fine)
            for (int i = 0; i < n_coarse_total; ++i) {
               in[i] = mtrandom::gaussian() * sigma;
            }

            // 2-4. Forward FFT → multiply sqrt(PSD) → inverse FFT
            run_fft_pipeline(in, forward, out, backward, result,
                             sqrt_PSD_coarse, n_coarse_total);

            // 5. Store with proper scaling
            for (int j = 0; j < n_coarse_total; ++j) {
               const size_t index = static_cast<size_t>(j) + static_cast<size_t>(r) * static_cast<size_t>(n_coarse_total);
               noise_field[index] = result[j] * norm_factor * inv_sqrt_S0 * scale;
            }

            // Progress bar (every 5%)
            int current_percent = static_cast<int>((r + 1) * 100.0 / realizations);
            if (current_percent >= last_printed_percent + 5 || r == realizations - 1) {
               float progress = static_cast<float>(r + 1) / realizations;
               int pos = static_cast<int>(bar_width * progress);
               std::cout << "\r  [";
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
         std::cout << std::endl;

         // Cleanup temporary FFTW resources
         fftw_destroy_plan(forward);
         fftw_destroy_plan(backward);
         fftw_free(in);
         fftw_free(out);
         fftw_free(result);

         #else
         std::cerr << "Error: quantum thermostat requires FFTW. Recompile with -DFFT -lfftw3." << std::endl;
         err::vexit();
         #endif
      }

      //=====================================================================
      // Linear interpolation on the coarse noise grid
      //
      // Given a fine step index (possibly fractional for RK4 sub-steps),
      // compute the corresponding coarse grid position and linearly
      // interpolate between the two bracketing coarse samples.
      //=====================================================================
      double get_noise(const std::vector<double>& coarse_noise,
                       double fine_step_idx, int M, size_t atom_idx) {

         double coarse_idx_float = fine_step_idx / M;
         size_t j = static_cast<size_t>(coarse_idx_float);
         double frac = coarse_idx_float - j;

         const size_t index1 = j + atom_idx;
         const size_t index2 = j + atom_idx + 1;

         // Bounds check: if the integrator ever asks for a sample past the
         // valid region, that's an upstream bug (likely a missed call to
         // update_noise_if_needed() or a miscount of window_start_fine).
         // Fail visibly: assert in debug; one-shot stderr message + clamp
         // in release so the run still produces output but the user sees
         // the flag.
         assert(index2 < coarse_noise.size()
                && "quantum::get_noise: index past end of coarse_noise buffer");
         if (index2 >= coarse_noise.size()) {
            static bool warned = false;
            if (!warned) {
               std::cerr << "Warning: quantum::get_noise indexed past end of "
                         << "coarse_noise buffer (index=" << index2
                         << ", size=" << coarse_noise.size() << "). "
                         << "Clamping; further occurrences silenced." << std::endl;
               warned = true;
            }
            return coarse_noise[coarse_noise.size() - 1];
         }

         return coarse_noise[index1] * (1.0 - frac) + coarse_noise[index2] * frac;
      }

   } // end of internal namespace

   //========================================================================
   // Public API: get_field
   //
   // Returns the quantum noise field for a given atom and spatial component
   // (0=x, 1=y, 2=z) at the current fine-step index. One sample per RK4
   // step — the same value is reused across K1-K4 (LSF_RK4 / HO convention).
   //
   // In windowed mode, the noise_index is offset relative to the window
   // start so that the lookup into the coarse_noise_field is correct.
   //========================================================================
   double get_field(int atom, int component) {
      using namespace internal;

      size_t idx = 0;
      if      (component == 0) idx = atom_idx_x[atom];
      else if (component == 1) idx = atom_idx_y[atom];
      else                     idx = atom_idx_z[atom];

      // Offset relative to window start (non-windowed mode: window_start_fine == 0).
      // One sample per LLG step — no within-step sub-stage variation.
      const double local_fine_idx = noise_index - window_start_fine;

      return get_noise(coarse_noise_field, local_fine_idx, M_decimation, idx);
   }

   //========================================================================
   // Public API: increment_time
   //
   // Advances the noise time index by one fine step.
   // Called at the end of each LLG integration step.
   //========================================================================
   void increment_time() {
      internal::noise_index += 1;
   }

   namespace internal{

      //=====================================================================
      // Export / append actual simulation noise for atom 0 to file
      //
      // Writes the x, y, z noise components for atom 0 from the current
      // coarse_noise_field buffer. On the first call the file header is
      // written (truncate mode); on subsequent calls data is appended.
      // This allows the full noise time series — including window
      // boundaries — to be inspected when using windowed generation.
      //=====================================================================
      void export_noise_data(const std::string& filename) {
         #ifdef FFT

         const int nc = valid_n_coarse;
         const double dt_coarse = stored_dt_fine * stored_M;

         // First call: write header (truncate). Subsequent calls: append.
         std::ios_base::openmode mode = noise_export_header_written
            ? (std::ios_base::out | std::ios_base::app)
            : std::ios_base::out;

         std::ofstream outfile(filename, mode);
         if(!outfile.is_open()){
            std::cerr << "Error: Unable to open noise export file " << filename << std::endl;
            err::vexit();
         }

         // Write header on first call only
         if(!noise_export_header_written){
            outfile << "# Quantum noise export — actual simulation noise for atom 0\n";
            outfile << "# noise_type: " << noise_type
                    << " (0=classical, 1=quantum, 2=quantum-no-zero)\n";
            outfile << "# windowed_mode: " << (windowed_mode ? "true" : "false") << "\n";
            outfile << "# valid_n_coarse: " << nc << "\n";
            if (windowed_mode) {
               outfile << "# fft_n_coarse: " << window_n_coarse << "\n";
               outfile << "# noise_seg: " << noise_seg << " (overlap-save segment size)\n";
            }
            outfile << "# dt_coarse: " << dt_coarse << " s\n";
            outfile << "# T_scaled: " << stored_T << "\n";
            outfile << "# Column 1: time (s)\n";
            outfile << "# Column 2: noise_x (atom 0)\n";
            outfile << "# Column 3: noise_y (atom 0)\n";
            outfile << "# Column 4: noise_z (atom 0)\n";
            noise_export_header_written = true;
         }

         // Atom 0 offsets into coarse_noise_field
         const size_t ix = atom_idx_x[0];
         const size_t iy = atom_idx_y[0];
         const size_t iz = atom_idx_z[0];

         outfile << std::scientific << std::setprecision(12);
         for (int k = 0; k < nc; ++k) {
            // stored_dt_fine is mp::dt, the REDUCED step (dt_SI * gamma).
            // Convert back to seconds so the exported axis matches its header
            // and a derived PSD lands on a real frequency scale.
            const double t = (window_start_fine + static_cast<double>(k) * stored_M)
                             * stored_dt_fine / mp::gamma_SI;
            outfile << t
                    << " " << coarse_noise_field[ix + k]
                    << " " << coarse_noise_field[iy + k]
                    << " " << coarse_noise_field[iz + k]
                    << "\n";
         }

         if(outfile.fail()){
            std::cerr << "Error: Failed to write noise data to " << filename << std::endl;
            err::vexit();
         }
         outfile.close();

         #else
         std::cerr << "Error: noise export requires FFTW." << std::endl;
         err::vexit();
         #endif
      }

      //=====================================================================
      // Initialize persistent noise structures for windowed mode
      //
      // Allocates the FFTW plans and arrays that persist across multiple
      // calls to generate_noise_window(). Precomputes scaling constants
      // and the sqrt(PSD) filter for the window frequency grid.
      //=====================================================================
      void init_noise_structures(int n_coarse, int realizations,
                                 double dt_fine, int M, double T) {
         #ifdef FFT

         // Store parameters for window regeneration
         window_n_coarse = n_coarse;
         num_realizations = realizations;
         stored_dt_fine = dt_fine;
         stored_M = M;
         stored_T = T;
         window_start_fine = 0;

         const double dt_coarse = dt_fine * M;

         // Precompute scaling constants
         const double S0 = material_S0_array[0];
         noise_inv_sqrt_S0 = (S0 > 0.0) ? 1.0 / std::sqrt(S0) : 1.0;
         noise_norm_factor = 1.0 / n_coarse;
         noise_sigma = 1.0 / std::sqrt(dt_fine);
         noise_scale = std::sqrt(dt_fine / dt_coarse);

         // Allocate persistent FFTW arrays
         fft_in     = (double*)fftw_malloc(sizeof(double) * n_coarse);
         fft_out    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n_coarse / 2 + 1));
         fft_result = (double*)fftw_malloc(sizeof(double) * n_coarse);

         // Create persistent FFTW plans (FFTW_MEASURE for best performance)
         fft_forward  = fftw_plan_dft_r2c_1d(n_coarse, fft_in, fft_out, FFTW_MEASURE);
         fft_backward = fftw_plan_dft_c2r_1d(n_coarse, fft_out, fft_result, FFTW_MEASURE);

         // Precompute sqrt(PSD) for the window frequency grid
         const double df = 1.0 / (n_coarse * dt_coarse);
         sqrt_PSD_window.resize(n_coarse / 2 + 1);
         for (int i = 0; i <= n_coarse / 2; ++i) {
            const double omega = 2.0 * M_PI * i * df;
            sqrt_PSD_window[i] = std::sqrt(PSD(omega, T, 0));
         }

         // Overlap-save segment parameters and buffer allocation
         if (windowed_mode) {
            noise_seg = n_coarse / OVERLAP_SAVE_SEGMENTS;
            valid_n_coarse = OVERLAP_SAVE_VALID * noise_seg;

            // Allocate overlap buffer: stores last N_OVERLAP segments of white
            // noise per realization for use as the next window's overlap input.
            // For the very first window there is no "previous" — fill with fresh
            // Gaussian samples at the same noise_sigma so the spectral shaping
            // of the first window sees a fully-populated white-noise input.
            // (Otherwise the convolution with the zero overlap leaves the first
            //  ~one-segment-worth of valid output near zero.)
            constexpr int N_OVERLAP_INIT = OVERLAP_SAVE_SEGMENTS - OVERLAP_SAVE_VALID; // = 2
            const size_t overlap_size = static_cast<size_t>(num_realizations) * N_OVERLAP_INIT * noise_seg;
            prev_white_noise.resize(overlap_size);
            for (size_t i = 0; i < overlap_size; ++i) {
               prev_white_noise[i] = mtrandom::gaussian() * noise_sigma;
            }

            // Allocate noise buffer for the valid region only
            const size_t total_elements = static_cast<size_t>(realizations) * static_cast<size_t>(valid_n_coarse);
            coarse_noise_field.resize(total_elements);

            std::cout << "  Windowed noise structures initialized (overlap-save):" << std::endl;
            std::cout << "    FFT coarse steps: " << n_coarse << std::endl;
            std::cout << "    Segment size: " << noise_seg << std::endl;
            std::cout << "    Valid coarse steps: " << valid_n_coarse << std::endl;
            std::cout << "    Realizations: " << realizations << std::endl;
            std::cout << "    Noise buffer: " << total_elements
                      << " (" << (total_elements * sizeof(double)) / (1024 * 1024) << " MB)" << std::endl;
            std::cout << "    Overlap buffer: " << overlap_size
                      << " (" << (overlap_size * sizeof(double)) / (1024 * 1024) << " MB)" << std::endl;
         } else {
            valid_n_coarse = n_coarse;
            noise_seg = 0;
            const size_t total_elements = static_cast<size_t>(realizations) * static_cast<size_t>(n_coarse);
            coarse_noise_field.resize(total_elements);

            std::cout << "  Noise structures initialized:" << std::endl;
            std::cout << "    Coarse steps: " << n_coarse << std::endl;
            std::cout << "    Realizations: " << realizations << std::endl;
            std::cout << "    Buffer size: " << total_elements
                      << " (" << (total_elements * sizeof(double)) / (1024 * 1024) << " MB)" << std::endl;
         }

         #else
         std::cerr << "Error: quantum thermostat requires FFTW." << std::endl;
         err::vexit();
         #endif
      }

      //=====================================================================
      // Generate one window of noise using persistent FFT resources
      //
      // Loops over all realizations (3 * num_atoms), generating independent
      // colored noise samples on the coarse grid. Uses the pre-allocated
      // FFTW plans and precomputed sqrt(PSD) filter.
      //=====================================================================
      void generate_noise_window() {
         #ifdef FFT

         constexpr int N_OVERLAP = OVERLAP_SAVE_SEGMENTS - OVERLAP_SAVE_VALID;  // = 2

         const int nc    = window_n_coarse;     // FFT length (OVERLAP_SAVE_SEGMENTS segments)
         const int seg   = noise_seg;            // Segment size (nc / OVERLAP_SAVE_SEGMENTS)
         const int valid = valid_n_coarse;       // Valid samples per window (OVERLAP_SAVE_VALID segments)

         for (int r = 0; r < num_realizations; ++r) {

            // 1. Copy previous window's white noise tail into the first
            //    N_OVERLAP segments (initialised to zero for the first window).
            const size_t overlap_base = static_cast<size_t>(r) * N_OVERLAP * seg;
            for (int i = 0; i < N_OVERLAP * seg; ++i) {
               fft_in[i] = prev_white_noise[overlap_base + i];
            }

            // 2. Fill the remaining OVERLAP_SAVE_VALID segments with fresh white noise.
            for (int i = N_OVERLAP * seg; i < nc; ++i) {
               fft_in[i] = mtrandom::gaussian() * noise_sigma;
            }

            // 3. Save the last N_OVERLAP segments of white noise for the next window
            //    (done before FFT to ensure fft_in is unmodified).
            for (int i = 0; i < N_OVERLAP * seg; ++i) {
               prev_white_noise[overlap_base + i] = fft_in[OVERLAP_SAVE_VALID * seg + i];
            }

            // 4-6. Forward FFT → multiply sqrt(PSD) → inverse FFT
            run_fft_pipeline(fft_in, fft_forward, fft_out,
                             fft_backward, fft_result,
                             sqrt_PSD_window, nc);

            // 7. Store the valid region: skip the first segment (corrupted by
            //    overlap with the prev window) and the last segment (corrupted
            //    by FFT wrap-around). The kept OVERLAP_SAVE_VALID segments
            //    sit at fft_result[seg .. (1+OVERLAP_SAVE_VALID)*seg).
            const size_t base = static_cast<size_t>(r) * static_cast<size_t>(valid);
            for (int j = 0; j < valid; ++j) {
               coarse_noise_field[base + j] = fft_result[seg + j] * noise_norm_factor * noise_inv_sqrt_S0 * noise_scale;
            }
         }

         // Append this window's noise for atom 0 if export is enabled (rank 0 only)
         if(export_noise){
            #ifdef MPICF
            if(vmpi::my_rank == 0) export_noise_data(export_noise_filename);
            #else
            export_noise_data(export_noise_filename);
            #endif
         }

         #endif
      }

      //=====================================================================
      // Check if the noise window needs regeneration
      //
      // In windowed mode, regenerate the noise buffer when the current
      // noise_index has consumed all fine steps in the active window that
      // can be safely interpolated.
      //
      // get_noise() does linear interpolation between coarse samples j
      // and j+1 where j = floor(local_fine_idx / M). For j+1 to stay
      // inside [0, valid_n_coarse), we need
      //   local_fine_idx < M · (valid_n_coarse − 1) = valid_fine_steps − M
      // so the last safe local index in a window is valid_fine_steps − M − 1.
      //
      // On regen we set window_start_fine = noise_index so the current
      // step's local index is 0 in the new buffer. This costs M fine
      // steps per window (the last M of the old buffer are never read)
      // but keeps the overlap-save bookkeeping simple and the assertion
      // in get_noise() satisfied.
      //=====================================================================
      void update_noise_if_needed() {
         if (!windowed_mode) return;

         const int valid_fine_steps = valid_n_coarse * stored_M;
         if ((noise_index - window_start_fine) < valid_fine_steps - stored_M) return;

         window_start_fine = noise_index;
         generate_noise_window();
      }

      //=====================================================================
      // Clean up persistent FFTW resources
      //=====================================================================
      void cleanup_noise_structures() {
         #ifdef FFT
         if (fft_forward)  { fftw_destroy_plan(fft_forward);  fft_forward  = nullptr; }
         if (fft_backward) { fftw_destroy_plan(fft_backward); fft_backward = nullptr; }
         if (fft_in)       { fftw_free(fft_in);               fft_in       = nullptr; }
         if (fft_out)      { fftw_free(fft_out);              fft_out      = nullptr; }
         if (fft_result)   { fftw_free(fft_result);           fft_result   = nullptr; }
         sqrt_PSD_window.clear();
         coarse_noise_field.clear();
         prev_white_noise.clear();
         #endif
      }

   } // end of internal namespace

   //=========================================================================
   // Public shutdown hook — releases the FFTW plans and buffers held by the
   // windowed-noise path. Safe to call even if the module was disabled.
   //=========================================================================
   void cleanup() {
      internal::cleanup_noise_structures();
   }

} // end of quantum namespace
