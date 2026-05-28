// src/DecayFunctions.cpp
#include "DecayFunctions.h"
#include "Constants.h"
#include "MathUtils.h"
#include "RadialGrid.h"
#include <cmath>

static double channel_Dirac(double mt, double y,
                            double mu_D, double M_baryon, double M_meson,
                            double branching, double stat_factor)
{
    if (mt <= M_baryon) return 0.0;
    double q = std::sqrt(mt*mt - M_baryon*M_baryon);
    double m = mD_central;
    double pstar = PStar(m, M_baryon, M_meson);
    if (pstar <= 0.0) return 0.0;

    double k_min = km(q, m, M_baryon, M_meson);
    double k_max = kp(q, m, M_baryon, M_meson);
    if (k_min >= k_max || k_min < 0.0) return 0.0;

    double factor = m / (4.0 * M_PI * pstar);
    double common = 1.0 / (mt * std::cosh(y) *
                    std::sqrt(mt*mt * std::cosh(y)*std::cosh(y) - M_baryon*M_baryon));

    const auto& grid = RadialGrid::instance();
    const int n_r = grid.r().size();
    const double dr = grid.dr();
    const auto& cosh_hr = grid.cosh();
    const auto& sinh_hr = grid.sinh();
    const auto& r2 = grid.r2();

    double integral = 0.0;
    for (int i = 0; i < n_r; ++i) {
        auto integrand_k = [&](double k) {
            double Ek = std::sqrt(m*m + k*k);
            double arg_exp = Ek * cosh_hr[i] * std::cosh(y) / T;
            double x = k * sinh_hr[i] / T;
            double boltz = std::exp(mu_D / T) * std::exp(-arg_exp) * sinh_over_x(x);
            return k * boltz * r2[i];
        };
        double integral_k = Integrate1D_high(integrand_k, k_min, k_max);
        integral += integral_k * dr;
    }
    return branching * stat_factor * factor * common * integral;
}

static double channel_PS(double mt, double y,
                         double mu_D, double M_baryon, double M_meson,
                         double branching, double stat_factor)
{
    if (mt <= M_baryon) return 0.0;
    double q = std::sqrt(mt*mt - M_baryon*M_baryon);
    double thresh = M_baryon + M_meson;

    const auto& grid = RadialGrid::instance();
    const int n_r = grid.r().size();
    const double dr = grid.dr();
    const auto& cosh_hr = grid.cosh();
    const auto& sinh_hr = grid.sinh();
    const auto& r2 = grid.r2();

    auto integrand_m = [&](double m) {
        double pstar = PStar(m, M_baryon, M_meson);
        if (pstar <= 0.0) return 0.0;

        double k_min = km(q, m, M_baryon, M_meson);
        double k_max = kp(q, m, M_baryon, M_meson);
        if (k_min >= k_max || k_min < 0.0) return 0.0;

        double factor = m / (4.0 * M_PI * pstar);
        double common = 1.0 / (mt * std::cosh(y) *
                        std::sqrt(mt*mt * std::cosh(y)*std::cosh(y) - M_baryon*M_baryon));

        double integral_rk = 0.0;
        for (int i = 0; i < n_r; ++i) {
            auto integrand_k = [&](double k) {
                double Ek = std::sqrt(m*m + k*k);
                double arg_exp = Ek * cosh_hr[i] * std::cosh(y) / T;
                double x = k * sinh_hr[i] / T;
                double boltz = std::exp(mu_D / T) * std::exp(-arg_exp) * sinh_over_x(x);
                return k * boltz * r2[i];
            };
            double integral_k = Integrate1D_high(integrand_k, k_min, k_max);
            integral_rk += integral_k * dr;
        }

        double ps = PhaseShiftWeight(m, mD_central, Gamma_Delta, thresh, mD_BW_max, M_baryon, M_meson);
        return ps * branching * stat_factor * factor * common * integral_rk;
    };

    return Integrate1D_high(integrand_m, thresh, mD_BW_max);
}

static double channel_BW(double mt, double y,
                         double mu_D, double M_baryon, double M_meson,
                         double branching, double stat_factor)
{
    if (mt <= M_baryon) return 0.0;
    double q = std::sqrt(mt*mt - M_baryon*M_baryon);
    double thresh = M_baryon + M_meson;

    const auto& grid = RadialGrid::instance();
    const int n_r = grid.r().size();
    const double dr = grid.dr();
    const auto& cosh_hr = grid.cosh();
    const auto& sinh_hr = grid.sinh();
    const auto& r2 = grid.r2();

    auto integrand_m = [&](double m) {
        double pstar = PStar(m, M_baryon, M_meson);
        if (pstar <= 0.0) return 0.0;

        double k_min = km(q, m, M_baryon, M_meson);
        double k_max = kp(q, m, M_baryon, M_meson);
        if (k_min >= k_max || k_min < 0.0) return 0.0;

        double factor = m / (4.0 * M_PI * pstar);
        double common = 1.0 / (mt * std::cosh(y) *
                        std::sqrt(mt*mt * std::cosh(y)*std::cosh(y) - M_baryon*M_baryon));

        double integral_rk = 0.0;
        for (int i = 0; i < n_r; ++i) {
            auto integrand_k = [&](double k) {
                double Ek = std::sqrt(m*m + k*k);
                double arg_exp = Ek * cosh_hr[i] * std::cosh(y) / T;
                double x = k * sinh_hr[i] / T;
                double boltz = std::exp(mu_D / T) * std::exp(-arg_exp) * sinh_over_x(x);
                return k * boltz * r2[i];
            };
            double integral_k = Integrate1D_high(integrand_k, k_min, k_max);
            integral_rk += integral_k * dr;
        }

        double bw = BreitWigner(m, mD_central, Gamma_Delta, thresh, mD_BW_max, M_baryon, M_meson);
        return bw * branching * stat_factor * factor * common * integral_rk;
    };

    return Integrate1D_high(integrand_m, thresh, mD_BW_max);
}

// ==========================================================================
// Protons – transverse mass spectra
// ==========================================================================
double dN_dmt_proton_from_Delta_Dirac(double mt, double y) {
    if (mt <= mp) return 0.0;
    double stat = g_Delta / g_proton;

    double res = 0.0;
    res += channel_Dirac(mt, y, mu_Dpp, mp, mpi, 1.0, stat);
    res += channel_Dirac(mt, y, mu_Dp,  mp, mpi0, 2.0/3.0, stat);
    res += channel_Dirac(mt, y, mu_D0,  mp, mpi,  1.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_proton_from_Delta_BW(double mt, double y) {
    if (mt <= mp) return 0.0;
    double stat = g_Delta / g_proton;

    double res = 0.0;
    res += channel_BW(mt, y, mu_Dpp, mp, mpi,   1.0,     stat);
    res += channel_BW(mt, y, mu_Dp,  mp, mpi0,  2.0/3.0, stat);
    res += channel_BW(mt, y, mu_D0,  mp, mpi,   1.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_proton_from_Delta_PS(double mt, double y) {
    if (mt <= mp) return 0.0;
    double stat = g_Delta / g_proton;

    double res = 0.0;
    res += channel_PS(mt, y, mu_Dpp, mp, mpi,   1.0,     stat);
    res += channel_PS(mt, y, mu_Dp,  mp, mpi0,  2.0/3.0, stat);
    res += channel_PS(mt, y, mu_D0,  mp, mpi,   1.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

// ==========================================================================
// Neutrons – transverse mass spectra
// ==========================================================================
double dN_dmt_neutron_from_Delta_Dirac(double mt, double y) {
    if (mt <= mn) return 0.0;
    double stat = g_Delta / g_neutron;

    double res = 0.0;
    res += channel_Dirac(mt, y, mu_Dp,  mn, mpi,   1.0/3.0, stat);
    res += channel_Dirac(mt, y, mu_D0,  mn, mpi0,  2.0/3.0, stat);
    res += channel_Dirac(mt, y, mu_Dm,  mn, mpi,   1.0,     stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_neutron_from_Delta_BW(double mt, double y) {
    if (mt <= mn) return 0.0;
    double stat = g_Delta / g_neutron;

    double res = 0.0;
    res += channel_BW(mt, y, mu_Dp,  mn, mpi,   1.0/3.0, stat);
    res += channel_BW(mt, y, mu_D0,  mn, mpi0,  2.0/3.0, stat);
    res += channel_BW(mt, y, mu_Dm,  mn, mpi,   1.0,     stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_neutron_from_Delta_PS(double mt, double y) {
    if (mt <= mn) return 0.0;
    double stat = g_Delta / g_neutron;

    double res = 0.0;
    res += channel_PS(mt, y, mu_Dp,  mn, mpi,   1.0/3.0, stat);
    res += channel_PS(mt, y, mu_D0,  mn, mpi0,  2.0/3.0, stat);
    res += channel_PS(mt, y, mu_Dm,  mn, mpi,   1.0,     stat);
    return res > 0.0 ? res : 0.0;
}

// ==========================================================================
// Positive pions – transverse mass spectra
// ==========================================================================
double dN_dmt_piplus_from_Delta_Dirac(double mt, double y) {
    if (mt <= mpi) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_Dirac(mt, y, mu_Dpp, mp, mpi,   1.0, stat);
    res += channel_Dirac(mt, y, mu_Dp,  mn, mpi,   1.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_piplus_from_Delta_BW(double mt, double y) {
    if (mt <= mpi) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_BW(mt, y, mu_Dpp, mp, mpi,   1.0,     stat);
    res += channel_BW(mt, y, mu_Dp,  mn, mpi,   1.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_piplus_from_Delta_PS(double mt, double y) {
    if (mt <= mpi) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_PS(mt, y, mu_Dpp, mp, mpi,   1.0,     stat);
    res += channel_PS(mt, y, mu_Dp,  mn, mpi,   1.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

// ==========================================================================
// Negative pions – transverse mass spectra
// ==========================================================================
double dN_dmt_piminus_from_Delta_Dirac(double mt, double y) {
    if (mt <= mpi) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_Dirac(mt, y, mu_D0,  mp, mpi,   1.0/3.0, stat);
    res += channel_Dirac(mt, y, mu_Dm,  mn, mpi,   1.0,     stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_piminus_from_Delta_BW(double mt, double y) {
    if (mt <= mpi) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_BW(mt, y, mu_D0,  mp, mpi,   1.0/3.0, stat);
    res += channel_BW(mt, y, mu_Dm,  mn, mpi,   1.0,     stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_piminus_from_Delta_PS(double mt, double y) {
    if (mt <= mpi) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_PS(mt, y, mu_D0,  mp, mpi,   1.0/3.0, stat);
    res += channel_PS(mt, y, mu_Dm,  mn, mpi,   1.0,     stat);
    return res > 0.0 ? res : 0.0;
}

// ==========================================================================
// Neutral pions – transverse mass spectra
// ==========================================================================
double dN_dmt_pi0_from_Delta_Dirac(double mt, double y) {
    if (mt <= mpi0) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_Dirac(mt, y, mu_Dp,  mp, mpi0,  2.0/3.0, stat);
    res += channel_Dirac(mt, y, mu_D0,  mn, mpi0,  2.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_pi0_from_Delta_BW(double mt, double y) {
    if (mt <= mpi0) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_BW(mt, y, mu_Dp,  mp, mpi0,  2.0/3.0, stat);
    res += channel_BW(mt, y, mu_D0,  mn, mpi0,  2.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

double dN_dmt_pi0_from_Delta_PS(double mt, double y) {
    if (mt <= mpi0) return 0.0;
    double stat = g_Delta / g_pion;

    double res = 0.0;
    res += channel_PS(mt, y, mu_Dp,  mp, mpi0,  2.0/3.0, stat);
    res += channel_PS(mt, y, mu_D0,  mn, mpi0,  2.0/3.0, stat);
    return res > 0.0 ? res : 0.0;
}

// ==========================================================================
// Rapidity distributions from Delta decays
// ==========================================================================
double dN_dy_proton_from_Delta_Dirac(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mp) return 0.0;
        return dN_dmt_proton_from_Delta_Dirac(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mp, mp + 5000.0);
}

double dN_dy_proton_from_Delta_BW(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mp) return 0.0;
        return dN_dmt_proton_from_Delta_BW(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mp, mp + 5000.0);
}

double dN_dy_proton_from_Delta_PS(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mp) return 0.0;
        return dN_dmt_proton_from_Delta_PS(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mp, mp + 5000.0);
}

double dN_dy_neutron_from_Delta_Dirac(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mn) return 0.0;
        return dN_dmt_neutron_from_Delta_Dirac(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mn, mn + 5000.0);
}

double dN_dy_neutron_from_Delta_BW(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mn) return 0.0;
        return dN_dmt_neutron_from_Delta_BW(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mn, mn + 5000.0);
}

double dN_dy_neutron_from_Delta_PS(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mn) return 0.0;
        return dN_dmt_neutron_from_Delta_PS(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mn, mn + 5000.0);
}

double dN_dy_piplus_from_Delta_Dirac(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi) return 0.0;
        return dN_dmt_piplus_from_Delta_Dirac(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi, mpi + 5000.0);
}

double dN_dy_piplus_from_Delta_BW(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi) return 0.0;
        return dN_dmt_piplus_from_Delta_BW(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi, mpi + 5000.0);
}

double dN_dy_piplus_from_Delta_PS(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi) return 0.0;
        return dN_dmt_piplus_from_Delta_PS(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi, mpi + 5000.0);
}

double dN_dy_piminus_from_Delta_Dirac(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi) return 0.0;
        return dN_dmt_piminus_from_Delta_Dirac(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi, mpi + 5000.0);
}

double dN_dy_piminus_from_Delta_BW(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi) return 0.0;
        return dN_dmt_piminus_from_Delta_BW(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi, mpi + 5000.0);
}

double dN_dy_piminus_from_Delta_PS(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi) return 0.0;
        return dN_dmt_piminus_from_Delta_PS(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi, mpi + 5000.0);
}

double dN_dy_pi0_from_Delta_Dirac(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi0) return 0.0;
        return dN_dmt_pi0_from_Delta_Dirac(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi0, mpi0 + 5000.0);
}

double dN_dy_pi0_from_Delta_BW(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi0) return 0.0;
        return dN_dmt_pi0_from_Delta_BW(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi0, mpi0 + 5000.0);
}

double dN_dy_pi0_from_Delta_PS(double y) {
    auto integrand = [&](double mt) {
        if (mt <= mpi0) return 0.0;
        return dN_dmt_pi0_from_Delta_PS(mt, y) * mt * mt;
    };
    return Integrate1D_high(integrand, mpi0, mpi0 + 5000.0);
}
