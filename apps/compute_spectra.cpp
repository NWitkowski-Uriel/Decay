#include <functional>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

struct Args {
    std::set<std::string> particles;
    std::string model = "all";
    std::string output = "spectra.root";
    int n_mt = 1000;
    int n_y = 101;
    double y_min = -2.0;
    double y_max = 2.0;
    double mt_max = 3000.0;
} args;

void parse_args(int argc, char* argv[]) {
    static option long_options[] = {
        {"particles", required_argument, nullptr, 'p'},
        {"model", required_argument, nullptr, 'm'},
        {"output", required_argument, nullptr, 'o'},
        {"nmt", required_argument, nullptr, 'n'},
        {"ny", required_argument, nullptr, 'y'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "p:m:o:n:y:", long_options, nullptr)) != -1) {
        switch (c) {
            case 'p': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    if (!item.empty()) args.particles.insert(item);
                }
                break;
            }
            case 'm': args.model = optarg; break;
            case 'o': args.output = optarg; break;
            case 'n': args.n_mt = std::max(50, std::atoi(optarg)); break;
            case 'y': args.n_y = std::max(11, std::atoi(optarg)); break;
            default: break;
        }
    }

    if (args.particles.empty()) {
        args.particles = {"proton", "neutron", "piplus", "piminus", "pi0"};
    }
    if (args.model == "both") args.model = "all";
}

bool want_model(const std::string& model) {
    return args.model == "all" || args.model == model;
}

struct ParticleInfo {
    double mass;
    double mu;
    double g;
    std::function<double(double,double)> dirac_mt;
    std::function<double(double,double)> bw_mt;
    std::function<double(double,double)> ps_mt;
    std::function<double(double)> dirac_y;
    std::function<double(double)> bw_y;
    std::function<double(double)> ps_y;
};

std::map<std::string, ParticleInfo> build_particle_map() {
    std::map<std::string, ParticleInfo> m;
    m["proton"] = {mp, mu_p, g_proton,
        dN_dmt_proton_from_Delta_Dirac, dN_dmt_proton_from_Delta_BW, dN_dmt_proton_from_Delta_PS,
        dN_dy_proton_from_Delta_Dirac, dN_dy_proton_from_Delta_BW, dN_dy_proton_from_Delta_PS};
    m["neutron"] = {mn, mu_n, g_neutron,
        dN_dmt_neutron_from_Delta_Dirac, dN_dmt_neutron_from_Delta_BW, dN_dmt_neutron_from_Delta_PS,
        dN_dy_neutron_from_Delta_Dirac, dN_dy_neutron_from_Delta_BW, dN_dy_neutron_from_Delta_PS};
    m["piplus"] = {mpi, mu_piplus, g_pion,
        dN_dmt_piplus_from_Delta_Dirac, dN_dmt_piplus_from_Delta_BW, dN_dmt_piplus_from_Delta_PS,
        dN_dy_piplus_from_Delta_Dirac, dN_dy_piplus_from_Delta_BW, dN_dy_piplus_from_Delta_PS};
    m["piminus"] = {mpi, mu_piminus, g_pion,
        dN_dmt_piminus_from_Delta_Dirac, dN_dmt_piminus_from_Delta_BW, dN_dmt_piminus_from_Delta_PS,
        dN_dy_piminus_from_Delta_Dirac, dN_dy_piminus_from_Delta_BW, dN_dy_piminus_from_Delta_PS};
    m["pi0"] = {mpi0, mu_pi0, g_pion,
        dN_dmt_pi0_from_Delta_Dirac, dN_dmt_pi0_from_Delta_BW, dN_dmt_pi0_from_Delta_PS,
        dN_dy_pi0_from_Delta_Dirac, dN_dy_pi0_from_Delta_BW, dN_dy_pi0_from_Delta_PS};
    return m;
}

void fill_mt_hist(TH1D& hist, double mass, const std::function<double(double)>& fn) {
    for (int i = 1; i <= hist.GetNbinsX(); ++i) {
        const double mt = hist.GetBinCenter(i);
        hist.SetBinContent(i, mt > mass ? fn(mt) : 0.0);
    }
}

void fill_y_hist(TH1D& hist, const std::function<double(double)>& fn) {
    for (int i = 1; i <= hist.GetNbinsX(); ++i) {
        const double y = hist.GetBinCenter(i);
        hist.SetBinContent(i, fn(y));
    }
}

} // namespace

int main(int argc, char* argv[]) {
    parse_args(argc, argv);
    RadialGrid::instance().setNumPoints(200);
    const auto particles = build_particle_map();

    TFile out(args.output.c_str(), "RECREATE");

    for (const auto& name : args.particles) {
        auto it = particles.find(name);
        if (it == particles.end()) continue;
        const auto& p = it->second;

        TH1D h_prim_mt(("h_prim_" + name + "_mt").c_str(), "", args.n_mt, p.mass, args.mt_max);
        TH1D h_prim_y(("h_prim_" + name + "_y").c_str(), "", args.n_y, args.y_min, args.y_max);

        fill_mt_hist(h_prim_mt, p.mass, [&](double mt) { return dN_dmt_primordial(mt, 0.0, p.mu, p.mass, p.g); });
        fill_y_hist(h_prim_y, [&](double y) { return dN_dy_primordial_full(y, p.mu, p.mass, p.g); });
        h_prim_mt.Write();
        h_prim_y.Write();

        auto write_model = [&](const std::string& model_name,
                               const std::function<double(double,double)>& f_mt,
                               const std::function<double(double)>& f_y) {
            TH1D h_delta_mt(("h_delta_" + model_name + "_" + name + "_mt").c_str(), "", args.n_mt, p.mass, args.mt_max);
            TH1D h_delta_y(("h_delta_" + model_name + "_" + name + "_y").c_str(), "", args.n_y, args.y_min, args.y_max);
            TH1D h_total_mt(("h_total_" + model_name + "_" + name + "_mt").c_str(), "", args.n_mt, p.mass, args.mt_max);
            TH1D h_total_y(("h_total_" + model_name + "_" + name + "_y").c_str(), "", args.n_y, args.y_min, args.y_max);

            fill_mt_hist(h_delta_mt, p.mass, [&](double mt) { return f_mt(mt, 0.0); });
            fill_y_hist(h_delta_y, [&](double y) { return f_y(y); });

            h_total_mt.Add(&h_prim_mt);
            h_total_mt.Add(&h_delta_mt);
            h_total_y.Add(&h_prim_y);
            h_total_y.Add(&h_delta_y);

            h_delta_mt.Write();
            h_delta_y.Write();
            h_total_mt.Write();
            h_total_y.Write();
        };

        if (want_model("dirac")) write_model("dirac", p.dirac_mt, p.dirac_y);
        if (want_model("bw")) write_model("bw", p.bw_mt, p.bw_y);
        if (want_model("ps")) write_model("ps", p.ps_mt, p.ps_y);
    }

    out.Close();
    std::cout << "Saved spectra to " << args.output << std::endl;
    return 0;
}
