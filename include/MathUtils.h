// include/MathUtils.h
#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <cmath>
#include <functional>

// ============================================
// Kinematic functions
// ============================================

/**
 * Compute the momentum of particle 1 in the rest frame of a resonance.
 * @param Mr  Resonance mass
 * @param M1  Mass of particle 1
 * @param M2  Mass of particle 2
 * @return    Momentum in the resonance rest frame (if kinematically allowed, else 0)
 */
double PStar(double Mr, double M1, double M2);

/**
 * Lower limit of resonance momentum integration (k_min) for a given final‑state particle momentum q.
 * Corresponds to formula (8) in the documentation.
 * @param q   Momentum of the final particle (e.g., proton)
 * @param Mr  Resonance mass
 * @param M1  Mass of the final particle
 * @param M2  Mass of the other decay product
 * @return    k_min
 */
double km(double q, double Mr, double M1, double M2);

/**
 * Upper limit of resonance momentum integration (k_max) for a given final‑state particle momentum q.
 * @param q   Momentum of the final particle
 * @param Mr  Resonance mass
 * @param M1  Mass of the final particle
 * @param M2  Mass of the other decay product
 * @return    k_max
 */
double kp(double q, double Mr, double M1, double M2);

// ============================================
// Breit‑Wigner distribution
// ============================================

/**
 * Normalized Breit‑Wigner distribution.
 * The distribution is normalized to unity over the interval [m_min, m_max].
 * @param m      Mass value
 * @param m0     Central mass
 * @param Gamma  Width
 * @param m_min  Lower integration bound (kinematic threshold)
 * @param m_max  Upper integration bound (e.g., 1500 MeV)
 * @return       Value of the normalized distribution at m
 */
double BreitWigner(double m, double m0, double Gamma, double m_min, double m_max);

// ============================================
// Numerical integration routines
// ============================================

/**
 * High‑accuracy 1D integration using the trapezoidal rule.
 * @param func  Function of one variable
 * @param a     Lower limit
 * @param b     Upper limit
 * @return      Approximate integral
 */
double Integrate1D_high(const std::function<double(double)>& func, double a, double b);

/**
 * High‑accuracy 2D integration using the trapezoidal rule (midpoint rule in each dimension).
 * @param func   Function of two variables (x,y)
 * @param x_min  Lower x limit
 * @param x_max  Upper x limit
 * @param y_min  Lower y limit
 * @param y_max  Upper y limit
 * @return       Approximate double integral
 */
double Integrate2D_high(const std::function<double(double,double)>& func,
                        double x_min, double x_max,
                        double y_min, double y_max);

// ============================================
// Auxiliary math helpers
// ============================================

/**
 * Safe computation of sinh(x)/x for small x (returns 1 for x ≈ 0).
 */
inline double sinh_over_x(double x) {
    const double eps = 1e-10;
    if (std::abs(x) < eps) return 1.0;
    return std::sinh(x) / x;
}

#endif // MATHUTILS_H