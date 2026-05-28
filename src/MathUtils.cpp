// src/MathUtils.cpp
#include "MathUtils.h"
#include "Constants.h"
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

namespace {

struct CacheKey {
    double m_min;
    double m_max;
    double m1;
    double m2;

    bool operator<(const CacheKey& other) const {
        return std::tie(m_min, m_max, m1, m2) < std::tie(other.m_min, other.m_max, other.m1, other.m2);
    }
};

struct SpectralTable {
    double x_min = 0.0;
    double x_max = 0.0;
    double dx = 0.0;
    std::vector<double> y;
};

constexpr int kSpectralTableSize = 2048;

static std::map<CacheKey, SpectralTable> bw_tables;
static std::map<CacheKey, SpectralTable> ps_tables;

double PStarPrime(double M, double m1, double m2) {
    const double ps = PStar(M, m1, m2);
    if (ps <= 0.0 || M <= 0.0) return 0.0;
    const double num = M*M*M*M - std::pow(m1*m1 - m2*m2, 2);
    return num / (4.0 * M*M*M * ps);
}

double bw_raw(double M) {
    return 1.0 / ((M - mD_central) * (M - mD_central) + 0.25 * Gamma_Delta * Gamma_Delta);
}

double PiPS(double M, double m1, double m2) {
    const double ps = PStar(M, m1, m2);
    return c1PS * ps * ps + c2PS * ps * ps * ps * ps;
}

double deltaArgPS(double M, double m1, double m2) {
    const double ps = PStar(M, m1, m2);
    if (ps <= 0.0 || M <= 0.0 || std::abs(M*M - mD_central*mD_central) < 1e-12) return 0.0;
    return -(2.0 * alpha0PS * ps * ps * ps) /
           (3.0 * M * (M*M - mD_central*mD_central) * (1.0 + PiPS(M, m1, m2)));
}

double dDeltaPS(double M, double m1, double m2) {
    const double ps = PStar(M, m1, m2);
    if (ps <= 0.0 || M <= m1 + m2 || M <= 0.0) return 0.0;

    const double dps = PStarPrime(M, m1, m2);
    const double y = std::atan(deltaArgPS(M, m1, m2));
    const double arg = deltaArgPS(M, m1, m2);
    const double denom = 1.0 + arg * arg;
    if (denom <= 0.0) return 0.0;

    const double pi_ps = PiPS(M, m1, m2);
    const double dPi = dps * (2.0 * c1PS * ps + 4.0 * c2PS * ps * ps * ps);
    const double bracket = (3.0 * dps / ps)
        - (1.0 / M)
        - (2.0 * M / (M*M - mD_central*mD_central))
        - (dPi / (1.0 + pi_ps));

    (void)y;
    return arg * bracket / denom;
}

double BPS(double M, double m1, double m2) {
    return 2.0 * dDeltaPS(M, m1, m2);
}

SpectralTable build_table(const CacheKey& key, const std::function<double(double)>& raw) {
    SpectralTable table;
    table.x_min = key.m_min;
    table.x_max = key.m_max;
    table.dx = (key.m_max - key.m_min) / (kSpectralTableSize - 1);
    table.y.resize(kSpectralTableSize, 0.0);

    for (int i = 0; i < kSpectralTableSize; ++i) {
        const double x = table.x_min + i * table.dx;
        const double v = raw(x);
        table.y[i] = (std::isfinite(v) && v > 0.0) ? v : 0.0;
    }

    auto interp = [&](double x) {
        if (x < table.x_min || x > table.x_max) return 0.0;
        double pos = (x - table.x_min) / table.dx;
        int i = static_cast<int>(pos);
        if (i < 0) return table.y.front();
        if (i >= kSpectralTableSize - 1) return table.y.back();
        double t = pos - i;
        return (1.0 - t) * table.y[i] + t * table.y[i + 1];
    };

    const double integral = Integrate1D_high(interp, key.m_min, key.m_max);
    const double norm = (integral > 0.0 && std::isfinite(integral)) ? 1.0 / integral : 1.0;
    for (double& v : table.y) v *= norm;
    return table;
}

const SpectralTable& get_bw_table(const CacheKey& key) {
    auto it = bw_tables.find(key);
    if (it == bw_tables.end()) {
        it = bw_tables.emplace(key, build_table(key, [](double x) { return bw_raw(x); })).first;
    }
    return it->second;
}

const SpectralTable& get_ps_table(const CacheKey& key) {
    auto it = ps_tables.find(key);
    if (it == ps_tables.end()) {
        it = ps_tables.emplace(key, build_table(key, [&](double x) { return BPS(x, key.m1, key.m2); })).first;
    }
    return it->second;
}

double eval_table(double x, const SpectralTable& table) {
    if (x < table.x_min || x > table.x_max) return 0.0;
    const double pos = (x - table.x_min) / table.dx;
    const int i = static_cast<int>(pos);
    if (i < 0) return table.y.front();
    if (i >= static_cast<int>(table.y.size()) - 1) return table.y.back();
    const double t = pos - i;
    const double value = (1.0 - t) * table.y[i] + t * table.y[i + 1];
    return (std::isfinite(value) && value > 0.0) ? value : 0.0;
}

} // namespace

double PStar(double Mr, double M1, double M2) {
    const double term1 = Mr*Mr - (M1 + M2)*(M1 + M2);
    const double term2 = Mr*Mr - (M1 - M2)*(M1 - M2);
    if (term1 <= 0.0 || term2 <= 0.0 || Mr <= 0.0) return 0.0;
    return std::sqrt(term1 * term2) / (2.0 * Mr);
}

double km(double q, double Mr, double M1, double M2) {
    const double pstar = PStar(Mr, M1, M2);
    if (pstar <= 1e-12) return 0.0;
    const double Estar = std::sqrt(M1*M1 + pstar*pstar);
    const double Eq    = std::sqrt(M1*M1 + q*q);
    return Mr * std::fabs(Estar*q - pstar*Eq) / (M1*M1);
}

double kp(double q, double Mr, double M1, double M2) {
    const double pstar = PStar(Mr, M1, M2);
    if (pstar <= 1e-12) return 0.0;
    const double Estar = std::sqrt(M1*M1 + pstar*pstar);
    const double Eq    = std::sqrt(M1*M1 + q*q);
    return Mr * std::fabs(Estar*q + pstar*Eq) / (M1*M1);
}

double BreitWigner(double M, double m0, double Gamma, double m_min, double m_max,
                   double m1, double m2) {
    (void)m0;
    (void)Gamma;
    const CacheKey key{m_min, m_max, m1, m2};
    return eval_table(M, get_bw_table(key));
}

double PhaseShiftWeight(double M, double m0, double Gamma, double m_min, double m_max,
                        double m1, double m2) {
    (void)m0;
    (void)Gamma;
    const CacheKey key{m_min, m_max, m1, m2};
    return eval_table(M, get_ps_table(key));
}

double Integrate1D_high(const std::function<double(double)>& func, double a, double b) {
    if (!(b > a)) return 0.0;
    const int n = 4000;
    const double h = (b - a) / n;
    double sum = 0.5 * (func(a) + func(b));
    for (int i = 1; i < n; ++i) sum += func(a + i * h);
    return sum * h;
}

double Integrate2D_high(const std::function<double(double,double)>& func,
                        double x_min, double x_max,
                        double y_min, double y_max) {
    if (!(x_max > x_min) || !(y_max > y_min)) return 0.0;
    const int nx = 160;
    const int ny = 160;
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
