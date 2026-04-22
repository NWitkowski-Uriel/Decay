// include/Thermodynamics.h
#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include "Constants.h"
#include "MathUtils.h"

// ============================================
// Primordial thermal distributions
// ============================================

double dN_dmt_primordial(double mt, double y, double mu, double M, double g = 1.0);
double dN_dy_primordial_full(double y, double mu, double M, double g = 1.0);

// ============================================
// Delta resonance primordial distributions
// ============================================

inline double dN_dmt_Delta_Dirac(double mt, double y, double mu, double g = 4.0) {
    return dN_dmt_primordial(mt, y, mu, mD_central, g);
}

inline double dN_dy_Delta_Dirac(double y, double mu, double g = 4.0) {
    return dN_dy_primordial_full(y, mu, mD_central, g);
}

double dN_dmt_Delta_BW(double mt, double y, double mu, double g = 4.0);
double dN_dy_Delta_BW(double y, double mu, double g = 4.0);

double dN_dmt_Delta_PS(double mt, double y, double mu, double g = 4.0);
double dN_dy_Delta_PS(double y, double mu, double g = 4.0);

#endif // THERMODYNAMICS_H
