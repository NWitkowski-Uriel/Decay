#include <cmath>
#include <fstream>
#include <filesystem>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "TFile.h"
#include "TTree.h"
#include "TColor.h"

#include "Constants.h"
#include "DecayFunctions.h"
#include "Plotting.h"
#include "RadialGrid.h"

namespace {

struct ParticleInfo {
    double mass;
    double mu;
    double g;
    std::function<double(double,double)> d_dirac_mt;
    std::function<double(double,double)> d_bw_mt;
    std::function<double(double,double)> d_ps_mt;
    std::function<double(double)> d_dirac_y;
    std::function<double(double)> d_bw_y;
    std::function<double(double)> d_ps_y;
};

std::map<std::string, ParticleInfo> particle_map() {
    return {
        {"proton", {mp, mu_p, g_proton, dN_dmt_proton_from_Delta_Dirac, dN_dmt_proton_from_Delta_BW, dN_dmt_proton_from_Delta_PS, dN_dy_proton_from_Delta_Dirac, dN_dy_proton_from_Delta_BW, dN_dy_proton_from_Delta_PS}},
        {"neutron", {mn, mu_n, g_neutron, dN_dmt_neutron_from_Delta_Dirac, dN_dmt_neutron_from_Delta_BW, dN_dmt_neutron_from_Delta_PS, dN_dy_neutron_from_Delta_Dirac, dN_dy_neutron_from_Delta_BW, dN_dy_neutron_from_Delta_PS}},
        {"piplus", {mpi, mu_piplus, g_pion, dN_dmt_piplus_from_Delta_Dirac, dN_dmt_piplus_from_Delta_BW, dN_dmt_piplus_from_Delta_PS, dN_dy_piplus_from_Delta_Dirac, dN_dy_piplus_from_Delta_BW, dN_dy_piplus_from_Delta_PS}},
        {"piminus", {mpi, mu_piminus, g_pion, dN_dmt_piminus_from_Delta_Dirac, dN_dmt_piminus_from_Delta_BW, dN_dmt_piminus_from_Delta_PS, dN_dy_piminus_from_Delta_Dirac, dN_dy_piminus_from_Delta_BW, dN_dy_piminus_from_Delta_PS}},
        {"pi0", {mpi0, mu_pi0, g_pion, dN_dmt_pi0_from_Delta_Dirac, dN_dmt_pi0_from_Delta_BW, dN_dmt_pi0_from_Delta_PS, dN_dy_pi0_from_Delta_Dirac, dN_dy_pi0_from_Delta_BW, dN_dy_pi0_from_Delta_PS}},
    };
}

}

int main() {
    RadialGrid::instance().setNumPoints(200);
    auto particles = particle_map();

    std::filesystem::create_directories("output");
    TFile out("output/step4c_analysis.root", "RECREATE");
    TTree table("multiplicity", "Integrated multiplicities by model");

    char particle[16];
    double primordial = 0, dirac = 0, bw = 0, ps = 0;
    table.Branch("particle", particle, "particle/C");
    table.Branch("primordial", &primordial);
    table.Branch("dirac_total", &dirac);
    table.Branch("bw_total", &bw);
    table.Branch("ps_total", &ps);

    std::ofstream csv("output/multiplicity_table.csv");
    csv << "particle,primordial,dirac_total,bw_total,ps_total\n";

    for (const auto& [name, p] : particles) {
        std::snprintf(particle, sizeof(particle), "%s", name.c_str());
        primordial = dN_dy_primordial_full(0.0, p.mu, p.mass, p.g);
        dirac = primordial + p.d_dirac_y(0.0);
        bw = primordial + p.d_bw_y(0.0);
        ps = primordial + p.d_ps_y(0.0);
        table.Fill();
        csv << name << ',' << primordial << ',' << dirac << ',' << bw << ',' << ps << '\n';

        std::vector<double> x(300), y_prim(300), y_d(300), y_b(300), y_p(300);
        const double mt0 = p.mass + 1e-3;
        const double step = 1000.0 / 299.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < 300; ++i) {
            const double mt = mt0 + i * step;
            x[i] = mt - p.mass;
            y_prim[i] = dN_dmt_primordial(mt, 0.0, p.mu, p.mass, p.g);
            y_d[i] = y_prim[i] + p.d_dirac_mt(mt, 0.0);
            y_b[i] = y_prim[i] + p.d_bw_mt(mt, 0.0);
            y_p[i] = y_prim[i] + p.d_ps_mt(mt, 0.0);
        }

        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> series = {
            {&x, &y_prim, kBlack, 1, "Primordial"},
            {&x, &y_d, kRed + 1, 1, "Dirac total"},
            {&x, &y_b, kBlue + 1, 1, "BW total"},
            {&x, &y_p, kGreen + 2, 1, "PS total"}
        };
        DrawComparisonPlot("output/step4c_" + name + "_mt", name + " mT", "m_{T}-m [MeV]", "1/m_{T}^{2} dN/(dm_{T} dy)", series, {}, true, 0, 1000, 1e-12, 1e-2);
    }

    table.Write();
    out.Close();
    std::cout << "Saved ROOT tables and plots in output/." << std::endl;
    return 0;
}
