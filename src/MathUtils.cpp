// src/MathUtils.cpp
#include "MathUtils.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct SpectralCacheKey {
    double m0;
    double gamma0;
    double mMin;
    double mMax;
    double m1;
    double m2;

    bool operator<(const SpectralCacheKey& other) const {
        return std::tie(m0, gamma0, mMin, mMax, m1, m2) <
               std::tie(other.m0, other.gamma0, other.mMin, other.mMax, other.m1, other.m2);
    }
};

double safe_pstar(double m, double m1, double m2) {
    return PStar(m, m1, m2);
}

double gamma_delta_running(double m, double m0, double gamma0, double m1, double m2) {
    if (m <= 0.0) return 0.0;
    const double p = safe_pstar(m, m1, m2);
    const double p0 = safe_pstar(m0, m1, m2);
    if (p <= 0.0 || p0 <= 0.0) return 0.0;

    const double ratio = p / p0;
    return gamma0 * (m0 / m) * ratio * ratio * ratio;
}

double breit_wigner_raw(double m, double m0, double gamma0, double m1, double m2) {
    const double gamma_m = gamma_delta_running(m, m0, gamma0, m1, m2);
    if (gamma_m <= 0.0) return 0.0;

    const double dm2 = m * m - m0 * m0;
    const double denom = dm2 * dm2 + (m0 * gamma_m) * (m0 * gamma_m);
    if (denom <= 0.0) return 0.0;

    const double value = (2.0 * m * m0 * gamma_m) / (kPi * denom);
    return std::isfinite(value) && value > 0.0 ? value : 0.0;
}

double phase_shift_delta(double m, double m0, double gamma0, double m1, double m2) {
    const double gamma_m = gamma_delta_running(m, m0, gamma0, m1, m2);
    if (gamma_m <= 0.0) return 0.0;

    const double real_part = m0 * m0 - m * m;
    const double imag_part = m0 * gamma_m;
    double delta = std::atan2(imag_part, real_part);
    if (delta < 0.0) delta += kPi;
    return delta;
}

double phase_shift_density_raw(double m, double m0, double gamma0, double mMin, double mMax,
                               double m1, double m2) {
    if (m < mMin || m > mMax) return 0.0;
    const double h = std::max(1e-3, 1e-4 * (mMax - mMin));
    const double ml = std::max(mMin, m - h);
    const double mr = std::min(mMax, m + h);
    if (mr <= ml) return 0.0;

    const double dl = phase_shift_delta(ml, m0, gamma0, m1, m2);
    const double dr = phase_shift_delta(mr, m0, gamma0, m1, m2);
    const double rho = (dr - dl) / (kPi * (mr - ml));
    return std::isfinite(rho) && rho > 0.0 ? rho : 0.0;
}

double normalized_weight(double m,
                         const SpectralCacheKey& key,
                         const std::function<double(double)>& raw) {
    if (m < key.mMin || m > key.mMax) return 0.0;

    static std::map<SpectralCacheKey, double> normCache;
    auto it = normCache.find(key);
    if (it == normCache.end()) {
        const double integral = Integrate1D_high(raw, key.mMin, key.mMax);
        const double norm = (integral > 0.0 && std::isfinite(integral)) ? 1.0 / integral : 1.0;
        it = normCache.emplace(key, norm).first;
    }

    const double value = it->second * raw(m);
    return std::isfinite(value) && value > 0.0 ? value : 0.0;
}

} // namespace

double PStar(double Mr, double M1, double M2) {
    const double term1 = Mr * Mr - (M1 + M2) * (M1 + M2);
    const double term2 = Mr * Mr - (M1 - M2) * (M1 - M2);
    if (term1 <= 0.0 || term2 <= 0.0 || Mr <= 0.0) return 0.0;
    return std::sqrt(term1 * term2) / (2.0 * Mr);
}

double km(double q, double Mr, double M1, double M2) {
    const double pstar = PStar(Mr, M1, M2);
    if (pstar <= 1e-12) return 0.0;
    const double Estar = std::sqrt(M1 * M1 + pstar * pstar);
    const double Eq = std::sqrt(M1 * M1 + q * q);
    return Mr * std::fabs(Estar * q - pstar * Eq) / (M1 * M1);
}

double kp(double q, double Mr, double M1, double M2) {
    const double pstar = PStar(Mr, M1, M2);
    if (pstar <= 1e-12) return 0.0;
    const double Estar = std::sqrt(M1 * M1 + pstar * pstar);
    const double Eq = std::sqrt(M1 * M1 + q * q);
    return Mr * std::fabs(Estar * q + pstar * Eq) / (M1 * M1);
}

double BreitWigner(double m, double m0, double Gamma, double m_min, double m_max,
                   double m1, double m2) {
    SpectralCacheKey key{m0, Gamma, m_min, m_max, m1, m2};
    auto raw = [&](double mm) { return breit_wigner_raw(mm, m0, Gamma, m1, m2); };
    return normalized_weight(m, key, raw);
}

double PhaseShiftWeight(double m, double m0, double Gamma, double m_min, double m_max,
                        double m1, double m2) {
    SpectralCacheKey key{m0, Gamma, m_min, m_max, m1, m2};
    auto raw = [&](double mm) { return phase_shift_density_raw(mm, m0, Gamma, m_min, m_max, m1, m2); };
    return normalized_weight(m, key, raw);
}

double Integrate1D_high(const std::function<double(double)>& func, double a, double b) {
    if (!(b > a)) return 0.0;
    int n = 4000;
    if (n % 2 != 0) ++n;
    const double h = (b - a) / n;

    double sum = func(a) + func(b);
    for (int i = 1; i < n; ++i) {
        const double x = a + i * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * func(x);
    }
    return sum * h / 3.0;
}

double Integrate2D_high(const std::function<double(double,double)>& func,
                        double x_min, double x_max,
                        double y_min, double y_max) {
    if (!(x_max > x_min) || !(y_max > y_min)) return 0.0;
    const int nx = 240;
    const int ny = 240;
    const double dx = (x_max - x_min) / nx;
    const double dy = (y_max - y_min) / ny;
    double sum = 0.0;

    for (int i = 0; i < nx; ++i) {
        const double x = x_min + (i + 0.5) * dx;
        for (int j = 0; j < ny; ++j) {
            const double y = y_min + (j + 0.5) * dy;
            sum += func(x, y);
        }
    }
    return sum * dx * dy;
}
