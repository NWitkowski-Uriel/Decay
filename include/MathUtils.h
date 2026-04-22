// include/MathUtils.h
#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <cmath>
#include <functional>

// ============================================
// Kinematic functions
// ============================================

double PStar(double Mr, double M1, double M2);
double km(double q, double Mr, double M1, double M2);
double kp(double q, double Mr, double M1, double M2);

// ============================================
// Spectral weights for Delta resonance
// ============================================

/**
 * Relativistic Breit-Wigner spectral function with energy-dependent p-wave width.
 * The function is normalized to unity on [m_min, m_max].
 */
double BreitWigner(double m, double m0, double Gamma, double m_min, double m_max,
                   double m1, double m2);

/**
 * Phase-shift spectral density rho_PS(m) = (1/pi) d delta / dm
 * constructed from the same energy-dependent p-wave width and normalized
 * to unity on [m_min, m_max].
 */
double PhaseShiftWeight(double m, double m0, double Gamma, double m_min, double m_max,
                        double m1, double m2);

// ============================================
// Numerical integration routines
// ============================================

double Integrate1D_high(const std::function<double(double)>& func, double a, double b);

double Integrate2D_high(const std::function<double(double,double)>& func,
                        double x_min, double x_max,
                        double y_min, double y_max);

// ============================================
// Auxiliary math helpers
// ============================================

inline double sinh_over_x(double x) {
    const double eps = 1e-10;
    if (std::abs(x) < eps) return 1.0;
    return std::sinh(x) / x;
}

#endif // MATHUTILS_H
