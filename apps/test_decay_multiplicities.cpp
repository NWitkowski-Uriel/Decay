#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

struct TestResult {
    std::string label;
    double value;
    double expected;
};

double integrate_y(const std::function<double(double)>& f, double y_max = 4.0, int n = 400) {
    double dy = 2.0 * y_max / n;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        double y = -y_max + i * dy;
        double w = (i == 0 || i == n) ? 0.5 : 1.0;
        sum += w * f(y);
    }
    return sum * dy;
}

double channel_dirac_mt(double mt, double y, double mu_D, double M_baryon, double M_meson,
                        double branching, double stat_factor) {
    if (mt <= M_baryon) return 0.0;
    double q = std::sqrt(mt * mt - M_baryon * M_baryon);
    double m = mD_central;
    double pstar = PStar(m, M_baryon, M_meson);
    if (pstar <= 0.0) return 0.0;

    double k_min = km(q, m, M_baryon, M_meson);
    double k_max = kp(q, m, M_baryon, M_meson);
    if (k_min >= k_max || k_min < 0.0) return 0.0;

    double factor = m / (4.0 * M_PI * pstar);
    double common = 1.0 / (mt * std::cosh(y) * std::sqrt(mt * mt * std::cosh(y) * std::cosh(y) - M_baryon * M_baryon));

    const auto& grid = RadialGrid::instance();
    double integral = 0.0;
    for (int i = 0; i < static_cast<int>(grid.r().size()); ++i) {
        auto integrand_k = [&](double k) {
            double Ek = std::sqrt(m * m + k * k);
            double arg_exp = Ek * grid.cosh()[i] * std::cosh(y) / T;
            double x = k * grid.sinh()[i] / T;
            double boltz = std::exp(mu_D / T) * std::exp(-arg_exp) * sinh_over_x(x);
            return k * boltz * grid.r2()[i];
        };
        integral += Integrate1D_high(integrand_k, k_min, k_max) * grid.dr();
    }
    return branching * stat_factor * factor * common * integral;
}

double channel_bw_mt(double mt, double y, double mu_D, double M_baryon, double M_meson,
                     double branching, double stat_factor) {
    if (mt <= M_baryon) return 0.0;
    double q = std::sqrt(mt * mt - M_baryon * M_baryon);
    double thresh = M_baryon + M_meson;

    const auto& grid = RadialGrid::instance();
    auto integrand_m = [&](double m) {
        double pstar = PStar(m, M_baryon, M_meson);
        if (pstar <= 0.0) return 0.0;
        double k_min = km(q, m, M_baryon, M_meson);
        double k_max = kp(q, m, M_baryon, M_meson);
        if (k_min >= k_max || k_min < 0.0) return 0.0;

        double factor = m / (4.0 * M_PI * pstar);
        double common = 1.0 / (mt * std::cosh(y) * std::sqrt(mt * mt * std::cosh(y) * std::cosh(y) - M_baryon * M_baryon));
        double integral_rk = 0.0;
        for (int i = 0; i < static_cast<int>(grid.r().size()); ++i) {
            auto integrand_k = [&](double k) {
                double Ek = std::sqrt(m * m + k * k);
                double arg_exp = Ek * grid.cosh()[i] * std::cosh(y) / T;
                double x = k * grid.sinh()[i] / T;
                double boltz = std::exp(mu_D / T) * std::exp(-arg_exp) * sinh_over_x(x);
                return k * boltz * grid.r2()[i];
            };
            integral_rk += Integrate1D_high(integrand_k, k_min, k_max) * grid.dr();
        }
        double bw = BreitWigner(m, mD_central, Gamma_Delta, thresh, mD_BW_max, M_baryon, M_meson);
        return bw * branching * stat_factor * factor * common * integral_rk;
    };
    return Integrate1D_high(integrand_m, thresh, mD_BW_max);
}

double channel_ps_mt(double mt, double y, double mu_D, double M_baryon, double M_meson,
                     double branching, double stat_factor) {
    if (mt <= M_baryon) return 0.0;
    double q = std::sqrt(mt * mt - M_baryon * M_baryon);
    double thresh = M_baryon + M_meson;

    const auto& grid = RadialGrid::instance();
    auto integrand_m = [&](double m) {
        double pstar = PStar(m, M_baryon, M_meson);
        if (pstar <= 0.0) return 0.0;
        double k_min = km(q, m, M_baryon, M_meson);
        double k_max = kp(q, m, M_baryon, M_meson);
        if (k_min >= k_max || k_min < 0.0) return 0.0;

        double factor = m / (4.0 * M_PI * pstar);
        double common = 1.0 / (mt * std::cosh(y) * std::sqrt(mt * mt * std::cosh(y) * std::cosh(y) - M_baryon * M_baryon));
        double integral_rk = 0.0;
        for (int i = 0; i < static_cast<int>(grid.r().size()); ++i) {
            auto integrand_k = [&](double k) {
                double Ek = std::sqrt(m * m + k * k);
                double arg_exp = Ek * grid.cosh()[i] * std::cosh(y) / T;
                double x = k * grid.sinh()[i] / T;
                double boltz = std::exp(mu_D / T) * std::exp(-arg_exp) * sinh_over_x(x);
                return k * boltz * grid.r2()[i];
            };
            integral_rk += Integrate1D_high(integrand_k, k_min, k_max) * grid.dr();
        }
        double ps = PhaseShiftWeight(m, mD_central, Gamma_Delta, thresh, mD_BW_max, M_baryon, M_meson);
        return ps * branching * stat_factor * factor * common * integral_rk;
    };
    return Integrate1D_high(integrand_m, thresh, mD_BW_max);
}

double channel_yield(const std::function<double(double,double,double,double,double,double,double)>& f,
                     double mu_D, double M_baryon, double M_meson, double branching, double stat_factor,
                     double mass_product) {
    return integrate_y([&](double y) {
        auto mt_integrand = [&](double mt) {
            return (mt > mass_product) ? f(mt, y, mu_D, M_baryon, M_meson, branching, stat_factor) * mt * mt : 0.0;
        };
        return Integrate1D_high(mt_integrand, mass_product, mass_product + 5000.0);
    });
}

std::string status(double value, double expected) {
    double rel = std::abs(value - expected) / (std::abs(expected) + 1e-15);
    if (!std::isfinite(value)) return "ERROR";
    if (rel < 1e-3) return "OK";
    if (rel < 1e-2) return "WARNING";
    return "ERROR";
}

void print_test_block(const std::string& title, const std::vector<TestResult>& tests) {
    std::cout << title << "\n";
    for (const auto& t : tests) {
        std::cout << "  " << std::left << std::setw(28) << t.label
                  << " = " << std::setw(12) << std::setprecision(8) << std::fixed << t.value
                  << " expected " << std::setw(12) << t.expected
                  << " [" << status(t.value, t.expected) << "]\n";
    }
}

void run_model(const std::string& model_name,
               const std::function<double(double,double,double,double,double,double,double)>& channel_mt,
               const std::function<double(double,double,double)>& delta_dy) {
    std::cout << "MODEL: " << model_name << "\n\n";

    const double stat_p = g_Delta / g_proton;
    const double stat_n = g_Delta / g_neutron;
    const double stat_pi = g_Delta / g_pion;

    double NDpp = integrate_y([&](double y){ return delta_dy(y, mu_Dpp, g_Delta); });
    double NDp  = integrate_y([&](double y){ return delta_dy(y, mu_Dp,  g_Delta); });
    double ND0  = integrate_y([&](double y){ return delta_dy(y, mu_D0,  g_Delta); });
    double NDm  = integrate_y([&](double y){ return delta_dy(y, mu_Dm,  g_Delta); });

    print_test_block("D++ -> p + pi+", {
        {"p / D++", channel_yield(channel_mt, mu_Dpp, mp, mpi, 1.0, stat_p, mp) / NDpp, 1.0},
        {"pi+ / D++", channel_yield(channel_mt, mu_Dpp, mp, mpi, 1.0, stat_pi, mpi) / NDpp, 1.0}
    });
    std::cout << "\n";

    print_test_block("D+ -> p + pi0, n + pi+", {
        {"p / ((2/3) D+)", channel_yield(channel_mt, mu_Dp, mp, mpi0, 2.0/3.0, stat_p, mp) / NDp, 2.0/3.0},
        {"pi0 / ((2/3) D+)", channel_yield(channel_mt, mu_Dp, mp, mpi0, 2.0/3.0, stat_pi, mpi0) / NDp, 2.0/3.0},
        {"n / ((1/3) D+)", channel_yield(channel_mt, mu_Dp, mn, mpi, 1.0/3.0, stat_n, mn) / NDp, 1.0/3.0},
        {"pi+ / ((1/3) D+)", channel_yield(channel_mt, mu_Dp, mn, mpi, 1.0/3.0, stat_pi, mpi) / NDp, 1.0/3.0}
    });
    std::cout << "\n";

    print_test_block("D0 -> n + pi0, p + pi-", {
        {"n / ((2/3) D0)", channel_yield(channel_mt, mu_D0, mn, mpi0, 2.0/3.0, stat_n, mn) / ND0, 2.0/3.0},
        {"pi0 / ((2/3) D0)", channel_yield(channel_mt, mu_D0, mn, mpi0, 2.0/3.0, stat_pi, mpi0) / ND0, 2.0/3.0},
        {"p / ((1/3) D0)", channel_yield(channel_mt, mu_D0, mp, mpi, 1.0/3.0, stat_p, mp) / ND0, 1.0/3.0},
        {"pi- / ((1/3) D0)", channel_yield(channel_mt, mu_D0, mp, mpi, 1.0/3.0, stat_pi, mpi) / ND0, 1.0/3.0}
    });
    std::cout << "\n";

    print_test_block("D- -> n + pi-", {
        {"n / D-", channel_yield(channel_mt, mu_Dm, mn, mpi, 1.0, stat_n, mn) / NDm, 1.0},
        {"pi- / D-", channel_yield(channel_mt, mu_Dm, mn, mpi, 1.0, stat_pi, mpi) / NDm, 1.0}
    });
    std::cout << "\n";
}

} // namespace

int main() {
    RadialGrid::instance().setNumPoints(200);
    std::cout << std::setprecision(8) << std::fixed;

    run_model("Dirac", channel_dirac_mt, dN_dy_Delta_Dirac);
    run_model("BW", channel_bw_mt, dN_dy_Delta_BW);
    run_model("PS", channel_ps_mt, dN_dy_Delta_PS);
    return 0;
}
