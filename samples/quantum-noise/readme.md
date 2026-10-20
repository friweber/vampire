# Quantum thermal noise

Atomistic spin dynamics of fcc Ni (3 nm cube, 2916 atoms) at 300 K with the
quantum noise module. The module replaces the classical white thermal field
by noise with a quantum spectral density (`quantum:noise-type`):

| File | Integrator | Noise |
|------|------------|-------|
| `input-classical` | llg-heun | classical white noise (reference) |
| `input` | llg-heun | quantum without zero point, `quantum-no-zero` |
| `input-quantum` | llg-heun | quantum with zero point, `quantum` |
| `input-open-system` | llg-quantum | open-system LLG: each spin coupled to a Lorentzian oscillator, `quantum-no-zero` |

```bash
cd samples/quantum-noise
../../vampire-serial                 # runs input
./run_simulations.sh                 # runs all four and prints the mean |m|
```

The quantum spectra suppress thermal noise above k_B T/hbar, so
`quantum-no-zero` gives a larger magnetisation than the classical run at the
same temperature. The zero-point term of `quantum` adds noise at all
temperatures and lowers the magnetisation of classical spins.

The open-system run needs the Lorentzian parameters in `Ni.mat`
(`quantum-lorentzian-central-frequency`, `quantum-lorentzian-width`, in
reduced units of 1.76e11 rad/s). Coloured noise also works with
`sim:integrator = spin-lattice`, where the lattice (force) noise is coloured
as well; see `src/quantum/readme` and the manual for all keywords.
