#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

using Clock = std::chrono::high_resolution_clock;

struct BenchCase {
    std::string name;
    std::function<double()> fn;
};

template <typename F>
double time_seconds(F&& f, int repeats) {
    volatile double sink = 0.0;
    const auto t0 = Clock::now();
    for (int i = 0; i < repeats; ++i) sink += f();
    const auto t1 = Clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    if (!std::isfinite(sink)) std::cerr << "Non-finite benchmark accumulator\n";
    return dt.count();
}

void run_case(const BenchCase& c, int repeats) {
    const double sec = time_seconds(c.fn, repeats);
    std::cout << std::left << std::setw(36) << c.name
              << "  total=" << std::setw(10) << std::fixed << std::setprecision(4) << sec
              << "  avg=" << std::setw(10) << (sec / repeats)
              << "  repeats=" << repeats << '\n';
}

} // namespace

int main() {
    RadialGrid::instance().setNumPoints(200);

#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << "\n";
#else
    std::cout << "OpenMP threads: disabled\n";
#endif

    const double y = 0.5;
    const double mt_p = mp + 500.0;
    const double mt_pi = mpi + 400.0;

    std::vector<BenchCase> cases = {
        {"primordial proton dN/dmt", [&]() { return dN_dmt_primordial(mt_p, y, mu_p, mp, g_proton); }},
        {"proton Delta Dirac dN/dmt", [&]() { return dN_dmt_proton_from_Delta_Dirac(mt_p, y); }},
        {"proton Delta BW dN/dmt", [&]() { return dN_dmt_proton_from_Delta_BW(mt_p, y); }},
        {"proton Delta PS dN/dmt", [&]() { return dN_dmt_proton_from_Delta_PS(mt_p, y); }},
        {"pion+ Delta BW dN/dmt", [&]() { return dN_dmt_piplus_from_Delta_BW(mt_pi, y); }},
        {"primordial proton dN/dy", [&]() { return dN_dy_primordial_full(y, mu_p, mp, g_proton); }},
        {"proton Delta BW dN/dy", [&]() { return dN_dy_proton_from_Delta_BW(y); }},
        {"proton Delta PS dN/dy", [&]() { return dN_dy_proton_from_Delta_PS(y); }},
    };

    std::cout << "\n=== Microbenchmarks ===\n";
    for (const auto& c : cases) {
        const int repeats = (c.name.find("dN/dy") != std::string::npos) ? 3 : 10;
        run_case(c, repeats);
    }

    std::cout << "\n=== Grid benchmark (export-like workload) ===\n";
    {
        const int n = 250;
        std::vector<double> mt(n);
        for (int i = 0; i < n; ++i) mt[i] = mp + 1e-3 + i * 4.0;

        const auto t0 = Clock::now();
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            volatile double a = dN_dmt_proton_from_Delta_BW(mt[i], 0.0);
            volatile double b = dN_dmt_proton_from_Delta_PS(mt[i], 0.0);
            (void)a; (void)b;
        }
        const auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "250-point BW+PS proton mt grid        total=" << std::fixed << std::setprecision(4) << dt.count() << " s\n";
    }

    return 0;
}
