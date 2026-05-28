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
const double mD_central = 1232.0;

// Intro.nb mass-distribution parameters
const double Gamma_Delta = 117.0;
const double mD_BW_min = mp + mpi;
const double mD_BW_max = 1500.0;

// Pok Man Lo et al. phase-shift parameters from Intro.nb
const double alpha0PS = 45.37;
const double c1PS = 16.7e-6;   // MeV^-2
const double c2PS = 65.6e-12;  // MeV^-4

// Spin degeneracy factors
const double g_Delta = 4.0;
const double g_proton = 2.0;
const double g_neutron = 2.0;
const double g_pion = 1.0;

// Chemical potentials for different particle species
const double mu_p      = muB + 0.5 * muI3;
const double mu_n      = muB - 0.5 * muI3;
const double mu_piplus = muI3;
const double mu_piminus = -muI3;
const double mu_pi0    = 0.0;

// Delta resonance chemical potentials for the four charge states
const double mu_Dpp = muB + 1.5 * muI3;
const double mu_Dp  = muB + 0.5 * muI3;
const double mu_D0  = muB - 0.5 * muI3;
const double mu_Dm  = muB - 1.5 * muI3;

// Experimental reference (for scaling)
const double N_proton_exp = 77.6;

#endif // CONSTANTS_H
