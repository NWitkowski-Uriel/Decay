// include/DecayFunctions.h
#ifndef DECAYFUNCTIONS_H
#define DECAYFUNCTIONS_H

#include "Constants.h"
#include "MathUtils.h"

// ============================================
// Decay contributions from Delta resonances
// (Differential spectra 1/mT^2 dN/(dmT dy))
// ============================================

// ---------- Protons ----------
/**
 * Proton spectrum from Delta decays (Dirac mass, fixed mD_central).
 * Includes all relevant charge states and branching ratios.
 */
double dN_dmt_proton_from_Delta_Dirac(double mt, double y);

/**
 * Proton spectrum from Delta decays with Breit‑Wigner mass smearing.
 * Integrates over resonance mass with proper kinematic thresholds.
 */
double dN_dmt_proton_from_Delta_BW(double mt, double y);

// ---------- Neutrons ----------
double dN_dmt_neutron_from_Delta_Dirac(double mt, double y);
double dN_dmt_neutron_from_Delta_BW(double mt, double y);

// ---------- Positive pions ----------
double dN_dmt_piplus_from_Delta_Dirac(double mt, double y);
double dN_dmt_piplus_from_Delta_BW(double mt, double y);

// ---------- Negative pions ----------
double dN_dmt_piminus_from_Delta_Dirac(double mt, double y);
double dN_dmt_piminus_from_Delta_BW(double mt, double y);

// ---------- Neutral pions ----------
double dN_dmt_pi0_from_Delta_Dirac(double mt, double y);
double dN_dmt_pi0_from_Delta_BW(double mt, double y);

// ============================================
// Rapidity distributions from Delta decays
// ============================================

// Protons
double dN_dy_proton_from_Delta_Dirac(double y);
double dN_dy_proton_from_Delta_BW(double y);

// Neutrons
double dN_dy_neutron_from_Delta_Dirac(double y);
double dN_dy_neutron_from_Delta_BW(double y);

// Positive pions
double dN_dy_piplus_from_Delta_Dirac(double y);
double dN_dy_piplus_from_Delta_BW(double y);

// Negative pions
double dN_dy_piminus_from_Delta_Dirac(double y);
double dN_dy_piminus_from_Delta_BW(double y);

// Neutral pions
double dN_dy_pi0_from_Delta_Dirac(double y);
double dN_dy_pi0_from_Delta_BW(double y);

#endif // DECAYFUNCTIONS_H