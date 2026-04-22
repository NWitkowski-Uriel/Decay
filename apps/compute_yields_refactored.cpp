#include <cmath>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

struct Args {
    std::set<std::string> particles;
    std::string model = "all";
    double y_max = 4.0;
    int n_y = 400;
} args;

void parse_args(int argc, char* argv[]) {
    static option long_options[] = {
        {"particles", required_argument, nullptr, 'p'},
        {"model", required_argument, nullptr, 'm'},
        {"ymax", required_argument, nullptr, 'y'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "p:m:y:", long_options, nullptr)) != -1) {
        switch (c) {
            case 'p': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ',')) args.particles.insert(item);
                break;
            }
            case 'm': args.model = optarg; break;
            case 'y': args.y_max = std::atof(optarg); break;
        }
    }
    if (args.particles.empty()) args.particles = {"proton"};
}

bool want(const std::string& m) {
    return args.model == "all" || args.model == m;
}

double integrate_y(const std::function<double(double)>& f) {
    double sum = 0.0;
    const double dy = 2.0 * args.y_max / args.n_y;
    for (int i = 0; i <= args.n_y; ++i) {
        double y = -args.y_max + i * dy;
        double w = (i == 0 || i == args.n_y) ? 0.5 : 1.0;
        sum += w * f(y);
    }
    return sum * dy;
}

} // namespace

int main(int argc, char* argv[]) {
    parse_args(argc, argv);
    RadialGrid::instance().setNumPoints(200);

    for (const auto& p : args.particles) {
        if (p == "proton") {
            double prim = integrate_y([&](double y){ return dN_dy_primordial_full(y, mu_p, mp, g_proton); });
            std::cout << "proton primordial = " << prim << std::endl;

            if (want("dirac")) {
                double d = integrate_y(dN_dy_proton_from_Delta_Dirac);
                std::cout << "dirac delta = " << d << " total = " << prim + d << std::endl;
            }
            if (want("bw")) {
                double d = integrate_y(dN_dy_proton_from_Delta_BW);
                std::cout << "bw delta = " << d << " total = " << prim + d << std::endl;
            }
            if (want("ps")) {
                double d = integrate_y(dN_dy_proton_from_Delta_PS);
                std::cout << "ps delta = " << d << " total = " << prim + d << std::endl;
            }

            if (!(std::isfinite(prim))) std::cerr << "ERROR: non-finite primordial" << std::endl;
        }
    }

    return 0;
}
