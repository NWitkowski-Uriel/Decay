// apps/compute_delta_decays.cpp
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <map>
#include <functional>
#include <sstream>
#include <set>

#include "TFile.h"
#include "TH1D.h"

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"

struct Args {
    std::string particles; // comma-separated list
    std::string model;     // "dirac", "bw", "ps", "all"
} args;

void parse(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"particles", required_argument, 0, 'p'},
        {"model",     required_argument, 0, 'm'},
        {0,0,0,0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "p:m:", long_options, nullptr)) != -1) {
        switch (c) {
            case 'p': args.particles = optarg; break;
            case 'm': args.model = optarg; break;
            default: break;
        }
    }

    if (args.particles.empty()) args.particles = "proton,neutron,piplus,piminus,pi0";
    if (args.model.empty()) args.model = "all";
}

std::set<std::string> parseParticleSet(const std::string& csv) {
    std::set<std::string> result;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) result.insert(item);
    }
    return result;
}

int main(int argc, char* argv[]) {
    parse(argc, argv);

    RadialGrid::instance().setNumPoints(200);

    const int n_mt = 1000;
    const double mt_min = 0.0;
    const double mt_max = 3000.0;

    const int n_y = 101;
    const double y_min = -2.0;
    const double y_max = 2.0;

    const auto selectedParticles = parseParticleSet(args.particles);
    const bool doDirac = (args.model == "dirac" || args.model == "all");
    const bool doBW = (args.model == "bw" || args.model == "all");
    const bool doPS = (args.model == "ps" || args.model == "all");

    if (!doDirac && !doBW && !doPS) {
        std::cerr << "Error: unsupported model '" << args.model << "'. Use dirac, bw, ps or all.\n";
        return 1;
    }

    std::map<std::string, double> particleMass = {
        {"proton", mp},
        {"neutron", mn},
        {"piplus", mpi},
        {"piminus", mpi},
        {"pi0", mpi0}
    };

    std::map<std::string, std::function<double(double,double)>> func_dirac_mt = {
        {"proton", dN_dmt_proton_from_Delta_Dirac},
        {"neutron", dN_dmt_neutron_from_Delta_Dirac},
        {"piplus", dN_dmt_piplus_from_Delta_Dirac},
        {"piminus", dN_dmt_piminus_from_Delta_Dirac},
        {"pi0", dN_dmt_pi0_from_Delta_Dirac}
    };

    std::map<std::string, std::function<double(double,double)>> func_bw_mt = {
        {"proton", dN_dmt_proton_from_Delta_BW},
        {"neutron", dN_dmt_neutron_from_Delta_BW},
        {"piplus", dN_dmt_piplus_from_Delta_BW},
        {"piminus", dN_dmt_piminus_from_Delta_BW},
        {"pi0", dN_dmt_pi0_from_Delta_BW}
    };

    std::map<std::string, std::function<double(double)>> func_dirac_y = {
        {"proton", dN_dy_proton_from_Delta_Dirac},
        {"neutron", dN_dy_neutron_from_Delta_Dirac},
        {"piplus", dN_dy_piplus_from_Delta_Dirac},
        {"piminus", dN_dy_piminus_from_Delta_Dirac},
        {"pi0", dN_dy_pi0_from_Delta_Dirac}
    };

    std::map<std::string, std::function<double(double)>> func_bw_y = {
        {"proton", dN_dy_proton_from_Delta_BW},
        {"neutron", dN_dy_neutron_from_Delta_BW},
        {"piplus", dN_dy_piplus_from_Delta_BW},
        {"piminus", dN_dy_piminus_from_Delta_BW},
        {"pi0", dN_dy_pi0_from_Delta_BW}
    };

    std::map<std::string, std::function<double(double,double)>> func_ps_mt = {
        {"proton", dN_dmt_proton_from_Delta_PS},
        {"neutron", dN_dmt_neutron_from_Delta_PS},
        {"piplus", dN_dmt_piplus_from_Delta_PS},
        {"piminus", dN_dmt_piminus_from_Delta_PS},
        {"pi0", dN_dmt_pi0_from_Delta_PS}
    };

    std::map<std::string, std::function<double(double)>> func_ps_y = {
        {"proton", dN_dy_proton_from_Delta_PS},
        {"neutron", dN_dy_neutron_from_Delta_PS},
        {"piplus", dN_dy_piplus_from_Delta_PS},
        {"piminus", dN_dy_piminus_from_Delta_PS},
        {"pi0", dN_dy_pi0_from_Delta_PS}
    };

    TFile f("delta.root", "RECREATE");
    std::vector<std::string> particles = {"proton", "neutron", "piplus", "piminus", "pi0"};

    for (const auto& p : particles) {
        if (!selectedParticles.empty() && selectedParticles.count(p) == 0) continue;

        const double m = particleMass[p];

        if (doDirac) {
            TH1D* h_mt = new TH1D(("h_delta_dirac_" + p + "_mt").c_str(), "", n_mt, mt_min, mt_max);
            TH1D* h_y  = new TH1D(("h_delta_dirac_" + p + "_y").c_str(), "", n_y, y_min, y_max);

            for (int i = 0; i < n_mt; ++i) {
                const double mt = mt_min + i * (mt_max - mt_min) / (n_mt - 1);
                if (mt <= m) continue;
                h_mt->SetBinContent(i + 1, func_dirac_mt[p](mt, 0.0));
            }
            for (int i = 0; i < n_y; ++i) {
                const double y = y_min + i * (y_max - y_min) / (n_y - 1);
                h_y->SetBinContent(i + 1, func_dirac_y[p](y));
            }

            h_mt->Write();
            h_y->Write();
            delete h_mt;
            delete h_y;
        }

        if (doBW) {
            TH1D* h_mt = new TH1D(("h_delta_bw_" + p + "_mt").c_str(), "", n_mt, mt_min, mt_max);
            TH1D* h_y  = new TH1D(("h_delta_bw_" + p + "_y").c_str(), "", n_y, y_min, y_max);

            for (int i = 0; i < n_mt; ++i) {
                const double mt = mt_min + i * (mt_max - mt_min) / (n_mt - 1);
                if (mt <= m) continue;
                h_mt->SetBinContent(i + 1, func_bw_mt[p](mt, 0.0));
            }
            for (int i = 0; i < n_y; ++i) {
                const double y = y_min + i * (y_max - y_min) / (n_y - 1);
                h_y->SetBinContent(i + 1, func_bw_y[p](y));
            }

            h_mt->Write();
            h_y->Write();
            delete h_mt;
            delete h_y;
        }

        if (doPS) {
            TH1D* h_mt = new TH1D(("h_delta_ps_" + p + "_mt").c_str(), "", n_mt, mt_min, mt_max);
            TH1D* h_y  = new TH1D(("h_delta_ps_" + p + "_y").c_str(), "", n_y, y_min, y_max);

            for (int i = 0; i < n_mt; ++i) {
                const double mt = mt_min + i * (mt_max - mt_min) / (n_mt - 1);
                if (mt <= m) continue;
                h_mt->SetBinContent(i + 1, func_ps_mt[p](mt, 0.0));
            }
            for (int i = 0; i < n_y; ++i) {
                const double y = y_min + i * (y_max - y_min) / (n_y - 1);
                h_y->SetBinContent(i + 1, func_ps_y[p](y));
            }

            h_mt->Write();
            h_y->Write();
            delete h_mt;
            delete h_y;
        }

        std::cout << "Computed Delta contributions for: " << p << "\n";
    }

    f.Close();
    return 0;
}
