// src/MathUtils.cpp
#include "MathUtils.h"
#include "RadialGrid.h"
#include <cmath>
#include <map>

// ============================================
// Kinematic functions (unchanged)
// ============================================
double PStar(double Mr, double M1, double M2) {
    double term1 = Mr*Mr - (M1 + M2)*(M1 + M2);
    double term2 = Mr*Mr - (M1 - M2)*(M1 - M2);
    if (term1 <= 0 || term2 <= 0) return 0.0;
    return std::sqrt(term1 * term2) / (2.0 * Mr);
}

double km(double q, double Mr, double M1, double M2) {
    double pstar = PStar(Mr, M1, M2);
    if (pstar <= 1e-10) return 0.0;
    double Estar = std::sqrt(M1*M1 + pstar*pstar);
    double Eq    = std::sqrt(M1*M1 + q*q);
    return Mr * std::fabs(Estar*q - pstar*Eq) / (M1*M1);
}

double kp(double q, double Mr, double M1, double M2) {
    double pstar = PStar(Mr, M1, M2);
    if (pstar <= 1e-10) return 0.0;
    double Estar = std::sqrt(M1*M1 + pstar*pstar);
    double Eq    = std::sqrt(M1*M1 + q*q);
    return Mr * std::fabs(Estar*q + pstar*Eq) / (M1*M1);
}

// ============================================
// Breit‑Wigner with caching
// ============================================
static double BreitWignerRaw(double m, double m0, double Gamma) {
    return (Gamma / (2.0 * M_PI)) / ((m - m0)*(m - m0) + (Gamma*0.5)*(Gamma*0.5));
}

double BreitWigner(double m, double m0, double Gamma, double m_min, double m_max) {
    if (m < m_min || m > m_max) return 0.0;
    // Use a map keyed by (m_min, m_max) to cache normalization
    static std::map<std::pair<double,double>, double> norm_cache;
    auto key = std::make_pair(m_min, m_max);
    auto it = norm_cache.find(key);
    if (it == norm_cache.end()) {
        auto integrand = [&](double mm) { return BreitWignerRaw(mm, m0, Gamma); };
        double I = Integrate1D_high(integrand, m_min, m_max);
        double norm = 1.0 / I;
        norm_cache[key] = norm;
        return norm * BreitWignerRaw(m, m0, Gamma);
    }
    return it->second * BreitWignerRaw(m, m0, Gamma);
}



static double PhaseShiftRaw(double m, double m0, double Gamma) {
    const double half_gamma = 0.5 * Gamma;
    const double x = (m - m0) / half_gamma;
    return (1.0 / M_PI) * std::atan(x) + 0.5;
}

double PhaseShiftWeight(double m, double m0, double Gamma, double m_min, double m_max) {
    if (m < m_min || m > m_max) return 0.0;
    static std::map<std::pair<double,double>, double> norm_cache;
    auto key = std::make_pair(m_min, m_max);
    auto it = norm_cache.find(key);
    if (it == norm_cache.end()) {
        auto integrand = [&](double mm) {
            const double half_gamma = 0.5 * Gamma;
            const double x = (mm - m0) / half_gamma;
            return 1.0 / (M_PI * half_gamma * (1.0 + x*x));
        };
        double I = Integrate1D_high(integrand, m_min, m_max);
        norm_cache[key] = (I > 0.0) ? (1.0 / I) : 1.0;
        return norm_cache[key] * integrand(m);
    }
    const double half_gamma = 0.5 * Gamma;
    const double x = (m - m0) / half_gamma;
    const double raw = 1.0 / (M_PI * half_gamma * (1.0 + x*x));
    return it->second * raw;
}
// ============================================
// Numerical integration – reduced versions
// ============================================
double Integrate1D_high(const std::function<double(double)>& func, double a, double b) {
    const int n = 2000;   // reduced for speed
    double h = (b - a) / n;
    double sum = 0.5 * (func(a) + func(b));
    for (int i = 1; i < n; ++i) sum += func(a + i * h);
    return sum * h;
}

double Integrate2D_high(const std::function<double(double,double)>& func,
                        double x_min, double x_max,
                        double y_min, double y_max) {
    const int nx = 80, ny = 80;
    double dx = (x_max - x_min) / nx;
    double dy = (y_max - y_min) / ny;
    double sum = 0.0;
    for (int i = 0; i < nx; ++i) {
        double x = x_min + (i + 0.5) * dx;
        for (int j = 0; j < ny; ++j) {
            double y = y_min + (j + 0.5) * dy;
            sum += func(x, y) * dx * dy;
        }
    }
    return sum;
}