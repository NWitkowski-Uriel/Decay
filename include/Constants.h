// include/Constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

// ============================================
// Physical constants and parameters
// ============================================

// Thermodynamic parameters
const double muB = 776.0;        // Baryon chemical potential [MeV]
const double muI3 = -14.1;       // Isospin chemical potential [MeV]
const double T   = 49.6;         // Temperature [MeV]
const double H   = 8.0;          // Hubble-like expansion rate [MeV]
const double R   = 16.02 / 197.0; // Fireball radius [MeV^-1] (converted from fm)

// Particle masses [MeV]
const double mp     = 938.27208816;   // Proton
const double mn     = 939.5654205;    // Neutron
const double mpi    = 139.57039;      // Charged pion
const double mpi0   = 134.9768;       // Neutral pion

// Delta resonance (central mass)
const double mD_central = 1232.0;      // PDG value [MeV]

// Breit‑Wigner parameters (used when finite width is considered)
const double Gamma_Delta = 120.0;       // Width of Delta [MeV]
const double mD_BW_min = mp + mpi;      // Lower integration limit for BW [MeV]
const double mD_BW_max   = 1500.0;      // Upper integration limit for BW [MeV]

// Spin degeneracy factors
const double g_Delta = 4.0;   // 2J+1 for J=3/2
const double g_proton = 2.0;  // 2J+1 for J=1/2
const double g_neutron = 2.0;
const double g_pion = 1.0;    // J=0

// Chemical potentials for different particle species
// (derived from muB and muI3 assuming isospin symmetry)
const double mu_p      = muB + 0.5 * muI3;      // Proton
const double mu_n      = muB - 0.5 * muI3;      // Neutron
const double mu_piplus = muI3;                   // π⁺
const double mu_piminus = -muI3;                  // π⁻
const double mu_pi0    = 0.0;                     // π⁰

// Delta resonance chemical potentials for the four charge states
const double mu_Dpp = muB + 1.5 * muI3;   // Δ⁺⁺
const double mu_Dp  = muB + 0.5 * muI3;   // Δ⁺
const double mu_D0  = muB - 0.5 * muI3;   // Δ⁰
const double mu_Dm  = muB - 1.5 * muI3;   // Δ⁻

// Experimental reference (for scaling)
const double N_proton_exp = 77.6;   // Number of free protons (used for scaling plots)

#endif // CONSTANTS_H