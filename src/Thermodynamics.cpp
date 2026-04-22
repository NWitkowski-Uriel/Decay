// src/Thermodynamics.cpp
#include "Thermodynamics.h"
#include "Constants.h"
#include "MathUtils.h"
#include "RadialGrid.h"
#include <cmath>

double dN_dmt_primordial(double mt, double y, double mu, double M, double g) {
    if (mt <= M) return 0.0;

    double E_p = mt * std::cosh(y);
    double pT  = std::sqrt(mt*mt - M*M);

    const auto& grid = RadialGrid::instance();
    const int n = grid.r().size();
    const double dr = grid.dr();
    const auto& cosh_hr = grid.cosh();
    const auto& sinh_hr = grid.sinh();
    const auto& r2 = grid.r2();

    double integral = 0.0;
    for (int i = 0; i < n; ++i) {
        double arg_exp = E_p * cosh_hr[i] / T;
        double x = pT * sinh_hr[i] / T;
        double boltz = std::exp(-arg_exp) * sinh_over_x(x);
        integral += r2[i] * boltz;
    }
    integral *= dr;

    double result = (g * std::exp(mu / T) * std::cosh(y) / (M_PI * pT)) * integral;
    return result > 0.0 ? result : 0.0;
}

double dN_dy_primordial_full(double y, double mu, double M, double g) {
    const auto& grid = RadialGrid::instance();
    const int n_r = grid.r().size();
    const double dr = grid.dr();
    const auto& cosh_hr = grid.cosh();
    const auto& sinh_hr = grid.sinh();
    const auto& r2 = grid.r2();

    double total = 0.0;
    for (int i = 0; i < n_r; ++i) {
        auto integrand_mt = [&](double mt) {
            if (mt <= M) return 0.0;
            double E_p = mt * std::cosh(y);
            double pT = std::sqrt(mt*mt - M*M);
            double arg_exp = E_p * cosh_hr[i] / T;
            double x = pT * sinh_hr[i] / T;
            double dNdmt_dy = (g * std::exp(mu / T) * std::cosh(y) / (M_PI * pT))
                              * r2[i] * std::exp(-arg_exp) * sinh_over_x(x);
            return dNdmt_dy * mt * mt;
        };
        double mt_min = M;
        double mt_max = M + 5000.0;
        double integral_mt = Integrate1D_high(integrand_mt, mt_min, mt_max);
        total += integral_mt * dr;
    }
    return total > 0.0 ? total : 0.0;
}

double dN_dmt_Delta_BW(double mt, double y, double mu, double g) {
    if (mt <= mD_BW_min) return 0.0;
    auto integrand = [&](double m) {
        if (mt <= m) return 0.0;
        double bw = BreitWigner(m, mD_central, Gamma_Delta, mD_BW_min, mD_BW_max, mp, mpi);
        double dnd = dN_dmt_primordial(mt, y, mu, m, g);
        return bw * dnd;
    };
    return Integrate1D_high(integrand, mD_BW_min, mD_BW_max);
}

double dN_dy_Delta_BW(double y, double mu, double g) {
    auto integrand_m = [&](double m) {
        auto integrand_mt_r = [&](double mt, double r) {
            if (mt <= m) return 0.0;
            double E_p = mt * std::cosh(y);
            double pT  = std::sqrt(mt*mt - m*m);
            double arg_exp = E_p * std::cosh(H * r) / T;
            double x = pT * std::sinh(H * r) / T;
            double dNdmt_dy = (g * std::exp(mu / T) * std::cosh(y) / (M_PI * pT))
                              * r * r * std::exp(-arg_exp) * sinh_over_x(x);
            return dNdmt_dy * mt * mt;
        };
        double res = Integrate2D_high(integrand_mt_r, m, m + 5000.0, 0.0, R);
        double bw = BreitWigner(m, mD_central, Gamma_Delta, mD_BW_min, mD_BW_max, mp, mpi);
        return bw * res;
    };
    return Integrate1D_high(integrand_m, mD_BW_min, mD_BW_max);
}

double dN_dmt_Delta_PS(double mt, double y, double mu, double g) {
    if (mt <= mD_BW_min) return 0.0;
    auto integrand = [&](double m) {
        if (mt <= m) return 0.0;
        const double ps = PhaseShiftWeight(m, mD_central, Gamma_Delta, mD_BW_min, mD_BW_max, mp, mpi);
        const double dnd = dN_dmt_primordial(mt, y, mu, m, g);
        return ps * dnd;
    };
    return Integrate1D_high(integrand, mD_BW_min, mD_BW_max);
}

double dN_dy_Delta_PS(double y, double mu, double g) {
    auto integrand_m = [&](double m) {
        auto integrand_mt_r = [&](double mt, double r) {
            if (mt <= m) return 0.0;
            double E_p = mt * std::cosh(y);
            double pT  = std::sqrt(mt*mt - m*m);
            double arg_exp = E_p * std::cosh(H * r) / T;
            double x = pT * std::sinh(H * r) / T;
            double dNdmt_dy = (g * std::exp(mu / T) * std::cosh(y) / (M_PI * pT))
                              * r * r * std::exp(-arg_exp) * sinh_over_x(x);
            return dNdmt_dy * mt * mt;
        };
        double res = Integrate2D_high(integrand_mt_r, m, m + 5000.0, 0.0, R);
        const double ps = PhaseShiftWeight(m, mD_central, Gamma_Delta, mD_BW_min, mD_BW_max, mp, mpi);
        return ps * res;
    };
    return Integrate1D_high(integrand_m, mD_BW_min, mD_BW_max);
}
