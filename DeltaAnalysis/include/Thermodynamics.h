// include/Thermodynamics.h
#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include "Constants.h"
#include "MathUtils.h"

// ============================================
// Primordial thermal distributions
// (Single freeze-out, Boltzmann statistics,
//  blast-wave with Hubble-like expansion)
// ============================================

/**
 * Compute the invariant transverse mass spectrum 1/mT^2 dN/(dmT dy)
 * for a primordial (directly emitted) particle.
 *
 * @param mt   Transverse mass (MeV)
 * @param y    Rapidity
 * @param mu   Chemical potential of the particle (MeV)
 * @param M    Rest mass of the particle (MeV)
 * @param g    Spin degeneracy factor (2J+1)
 * @return     1/mT^2 dN/(dmT dy) [MeV^{-3}]
 */
double dN_dmt_primordial(double mt, double y, double mu, double M, double g = 1.0);

/**
 * Compute the rapidity distribution dN/dy for a primordial particle,
 * obtained by integrating over transverse mass.
 *
 * @param y    Rapidity
 * @param mu   Chemical potential (MeV)
 * @param M    Rest mass (MeV)
 * @param g    Spin degeneracy
 * @return     dN/dy
 */
double dN_dy_primordial_full(double y, double mu, double M, double g = 1.0);

// ============================================
// Delta resonance primordial distributions
// (Dirac delta in mass: fixed mass mD_central)
// ============================================

/**
 * 1/mT^2 dN/(dmT dy) for Delta resonances with fixed mass (Dirac).
 * @param mt   Transverse mass
 * @param y    Rapidity
 * @param mu   Chemical potential of the specific Delta charge state
 * @param g    Spin degeneracy (usually 4)
 * @return     Spectrum
 */
inline double dN_dmt_Delta_Dirac(double mt, double y, double mu, double g = 4.0) {
    return dN_dmt_primordial(mt, y, mu, mD_central, g);
}

/**
 * dN/dy for Delta resonances with fixed mass (Dirac).
 */
inline double dN_dy_Delta_Dirac(double y, double mu, double g = 4.0) {
    return dN_dy_primordial_full(y, mu, mD_central, g);
}

// ============================================
// Delta resonance primordial distributions
// (Breit‑Wigner smearing over mass)
// ============================================

/**
 * 1/mT^2 dN/(dmT dy) for Delta resonances with mass distribution
 * given by a normalized Breit‑Wigner.
 *
 * @param mt   Transverse mass
 * @param y    Rapidity
 * @param mu   Chemical potential of the specific Delta charge state
 * @param g    Spin degeneracy (usually 4)
 * @return     Spectrum
 */
double dN_dmt_Delta_BW(double mt, double y, double mu, double g = 4.0);

/**
 * dN/dy for Delta resonances with Breit‑Wigner mass distribution.
 */
double dN_dy_Delta_BW(double y, double mu, double g = 4.0);

#endif // THERMODYNAMICS_H