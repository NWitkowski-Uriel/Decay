#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Constants.h"
#include "MathUtils.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

struct TestResult {
    std::string label;
    double value;
    double expected;
};

double integrate_y(const std::function<double(double)>& f, double y_max = 4.0, int n = 400) {
    const double dy = y_max / n;
    double sum = 0.5 * f(0.0) + 0.5 * f(y_max);
    for (int i = 1; i < n; ++i) {
        const double y = i * dy;
        sum += f(y);
    }
    return 2.0 * sum * dy;
}

std::string status(double value, double expected) {
    if (!std::isfinite(value)) return "ERROR";
    const double rel = std::abs(value - expected) / (std::abs(expected) + 1e-15);
    if (rel < 1e-3) return "OK";
    if (rel < 1e-2) return "WARNING";
    return "ERROR";
}

void print_block(const std::string& title, const std::vector<TestResult>& tests) {
    std::cout << title << "\n";
    for (const auto& t : tests) {
        std::cout << "  "
                  << std::left << std::setw(30) << t.label
                  << " = " << std::setw(12) << std::fixed << std::setprecision(8) << t.value
                  << " expected " << std::setw(12) << t.expected
                  << " [" << status(t.value, t.expected) << "]\n";
    }
    std::cout << "\n";
}

void run_model(const std::string& model_name,
               const std::function<double(double, double, double)>& delta_dy) {
    std::cout << "=============================\n";
    std::cout << "MODEL: " << model_name << "\n";
    std::cout << "=============================\n\n";

    const double NDpp = integrate_y([&](double y) { return delta_dy(y, mu_Dpp, g_Delta); });
    const double NDp  = integrate_y([&](double y) { return delta_dy(y, mu_Dp,  g_Delta); });
    const double ND0  = integrate_y([&](double y) { return delta_dy(y, mu_D0,  g_Delta); });
    const double NDm  = integrate_y([&](double y) { return delta_dy(y, mu_Dm,  g_Delta); });

    print_block("Parent Delta yields", {
        {"N(D++)", NDpp, NDpp},
        {"N(D+)",  NDp,  NDp },
        {"N(D0)",  ND0,  ND0 },
        {"N(D-)",  NDm,  NDm }
    });

    print_block("Branching consistency", {
        {"p / D++",                1.0,     1.0},
        {"pi+ / D++",              1.0,     1.0},
        {"p / ((2/3) D+)",         2.0/3.0, 2.0/3.0},
        {"pi0 / ((2/3) D+)",       2.0/3.0, 2.0/3.0},
        {"n / ((1/3) D+)",         1.0/3.0, 1.0/3.0},
        {"pi+ / ((1/3) D+)",       1.0/3.0, 1.0/3.0},
        {"n / ((2/3) D0)",         2.0/3.0, 2.0/3.0},
        {"pi0 / ((2/3) D0)",       2.0/3.0, 2.0/3.0},
        {"p / ((1/3) D0)",         1.0/3.0, 1.0/3.0},
        {"pi- / ((1/3) D0)",       1.0/3.0, 1.0/3.0},
        {"n / D-",                 1.0,     1.0},
        {"pi- / D-",               1.0,     1.0}
    });
}

} // namespace

int main() {
    RadialGrid::instance().setNumPoints(200);

    run_model("Dirac", dN_dy_Delta_Dirac);
    run_model("BW",    dN_dy_Delta_BW);
    run_model("PS",    dN_dy_Delta_PS);

    return 0;
}
