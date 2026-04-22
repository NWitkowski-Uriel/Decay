// apps/compute_primordial.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "TFile.h"
#include "TH1D.h"
#include "Constants.h"
#include "MathUtils.h"
#include "Thermodynamics.h"
#include "RadialGrid.h"

int main() {
    std::cout << "Initializing radial grid...\n";
    RadialGrid::instance().setNumPoints(200);
    std::cout << "Radial grid initialized with " << RadialGrid::instance().r().size() << " points.\n";

    // Ustaw mniejszą liczbę punktów dla szybszego testu (możesz później zwiększyć)
    const int n_mt = 3000;     // wcześniej 1000
    const int n_y = 101;       // wcześniej 101
    double mt_min = 0.0, mt_max = 3000.0;
    double y_min = -2.0, y_max = 2.0;

    std::cout << "Creating ROOT file primordial.root...\n";
    TFile f("primordial.root", "RECREATE");

    struct Particle {
        std::string name; double mass; double g; double mu;
    };
    std::vector<Particle> particles = {
        {"proton",  mp,  g_proton,  mu_p},
        {"neutron", mn,  g_neutron, mu_n},
        {"piplus",  mpi, g_pion,    mu_piplus},
        {"piminus", mpi, g_pion,    mu_piminus},
        {"pi0",     mpi0,g_pion,    mu_pi0}
    };

    for (const auto& p : particles) {
        std::cout << "\nProcessing particle: " << p.name << std::endl;

        // Histogram dla mt (y=0)
        TH1D *h_mt = new TH1D(("h_prim_" + p.name + "_mt").c_str(),
                               ("Primordial " + p.name + " m_T distribution").c_str(),
                               n_mt, mt_min, mt_max);
        // Histogram dla y
        TH1D *h_y = new TH1D(("h_prim_" + p.name + "_y").c_str(),
                             ("Primordial " + p.name + " rapidity distribution").c_str(),
                             n_y, y_min, y_max);

        // Obliczenia dla mt (y=0)
        std::cout << "  Computing mt distribution (y=0)...\n";
        #pragma omp parallel for
        for (int i = 0; i < n_mt; ++i) {
            double mt = mt_min + i * (mt_max - mt_min) / (n_mt - 1);
            if (mt <= p.mass) continue;
            double val = dN_dmt_primordial(mt, 0.0, p.mu, p.mass, p.g);
            h_mt->SetBinContent(i+1, val);
            if (i % (n_mt/10) == 0 && i != 0) {
                #pragma omp critical
                std::cout << "    mt progress: " << i << " / " << n_mt << std::endl;
            }
        }

        // Obliczenia dla y (całkowanie po mt)
        std::cout << "  Computing rapidity distribution...\n";
        #pragma omp parallel for
        for (int i = 0; i < n_y; ++i) {
            double y = y_min + i * (y_max - y_min) / (n_y - 1);
            double val = dN_dy_primordial_full(y, p.mu, p.mass, p.g);
            h_y->SetBinContent(i+1, val);
            if (i % (n_y/10) == 0 && i != 0) {
                #pragma omp critical
                std::cout << "    y progress: " << i << " / " << n_y << std::endl;
            }
        }

        h_mt->Write();
        h_y->Write();
        delete h_mt;
        delete h_y;
        std::cout << "  Finished " << p.name << std::endl;
    }

    f.Close();
    std::cout << "\nDone. Output saved to primordial.root\n";
    return 0;
}