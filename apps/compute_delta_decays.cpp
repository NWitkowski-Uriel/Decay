// compute_delta_decays.cpp
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <sstream>
#include <set>
#include <map>
#include <functional>

#include "TFile.h"
#include "TH1D.h"

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"

struct Args {
    std::set<std::string> particles;
    std::string model;     // dirac|bw|ps|all (both accepted as alias)
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
            case 'p': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    if (!item.empty()) args.particles.insert(item);
                }
                break;
            }
            case 'm': args.model = optarg; break;
        }
    }
    if (args.particles.empty()) args.particles = {"proton","neutron","piplus","piminus","pi0"};
    if (args.model.empty()) args.model = "all";
    if (args.model == "both") args.model = "all";
}

bool needParticle(const std::string& name) {
    return args.particles.count(name) > 0;
}

bool needModel(const std::string& m) {
    if (args.model == "all") return true;
    return args.model == m;
}

int main(int argc, char* argv[]) {
    parse(argc, argv);

    RadialGrid::instance().setNumPoints(200);

    const int n_mt = 1000;
    const double mt_min = 0.0, mt_max = 3000.0;
    const int n_y = 101;
    const double y_min = -2.0, y_max = 2.0;

    TFile f("delta.root", "RECREATE");

    std::map<std::string, std::function<double(double,double)>> func_dirac_mt;
    std::map<std::string, std::function<double(double,double)>> func_bw_mt;
    std::map<std::string, std::function<double(double,double)>> func_ps_mt;
    std::map<std::string, std::function<double(double)>> func_dirac_y;
    std::map<std::string, std::function<double(double)>> func_bw_y;
    std::map<std::string, std::function<double(double)>> func_ps_y;
    std::map<std::string, double> particle_mass;

    func_dirac_mt["proton"] = dN_dmt_proton_from_Delta_Dirac;
    func_bw_mt["proton"] = dN_dmt_proton_from_Delta_BW;
    func_ps_mt["proton"] = dN_dmt_proton_from_Delta_PS;
    func_dirac_y["proton"] = dN_dy_proton_from_Delta_Dirac;
    func_bw_y["proton"] = dN_dy_proton_from_Delta_BW;
    func_ps_y["proton"] = dN_dy_proton_from_Delta_PS;
    particle_mass["proton"] = mp;

    func_dirac_mt["neutron"] = dN_dmt_neutron_from_Delta_Dirac;
    func_bw_mt["neutron"] = dN_dmt_neutron_from_Delta_BW;
    func_ps_mt["neutron"] = dN_dmt_neutron_from_Delta_PS;
    func_dirac_y["neutron"] = dN_dy_neutron_from_Delta_Dirac;
    func_bw_y["neutron"] = dN_dy_neutron_from_Delta_BW;
    func_ps_y["neutron"] = dN_dy_neutron_from_Delta_PS;
    particle_mass["neutron"] = mn;

    func_dirac_mt["piplus"] = dN_dmt_piplus_from_Delta_Dirac;
    func_bw_mt["piplus"] = dN_dmt_piplus_from_Delta_BW;
    func_ps_mt["piplus"] = dN_dmt_piplus_from_Delta_PS;
    func_dirac_y["piplus"] = dN_dy_piplus_from_Delta_Dirac;
    func_bw_y["piplus"] = dN_dy_piplus_from_Delta_BW;
    func_ps_y["piplus"] = dN_dy_piplus_from_Delta_PS;
    particle_mass["piplus"] = mpi;

    func_dirac_mt["piminus"] = dN_dmt_piminus_from_Delta_Dirac;
    func_bw_mt["piminus"] = dN_dmt_piminus_from_Delta_BW;
    func_ps_mt["piminus"] = dN_dmt_piminus_from_Delta_PS;
    func_dirac_y["piminus"] = dN_dy_piminus_from_Delta_Dirac;
    func_bw_y["piminus"] = dN_dy_piminus_from_Delta_BW;
    func_ps_y["piminus"] = dN_dy_piminus_from_Delta_PS;
    particle_mass["piminus"] = mpi;

    func_dirac_mt["pi0"] = dN_dmt_pi0_from_Delta_Dirac;
    func_bw_mt["pi0"] = dN_dmt_pi0_from_Delta_BW;
    func_ps_mt["pi0"] = dN_dmt_pi0_from_Delta_PS;
    func_dirac_y["pi0"] = dN_dy_pi0_from_Delta_Dirac;
    func_bw_y["pi0"] = dN_dy_pi0_from_Delta_BW;
    func_ps_y["pi0"] = dN_dy_pi0_from_Delta_PS;
    particle_mass["pi0"] = mpi0;

    const std::vector<std::string> particles = {"proton","neutron","piplus","piminus","pi0"};
    for (const auto& p : particles) {
        if (!needParticle(p)) continue;

        auto write_model = [&](const std::string& model_name,
                               const std::function<double(double,double)>& f_mt,
                               const std::function<double(double)>& f_y) {
            TH1D *h_mt = new TH1D(("h_delta_" + model_name + "_" + p + "_mt").c_str(), "", n_mt, mt_min, mt_max);
            TH1D *h_y  = new TH1D(("h_delta_" + model_name + "_" + p + "_y").c_str(),  "", n_y, y_min, y_max);
            for (int i=0; i<n_mt; ++i) {
                double mt = mt_min + i*(mt_max-mt_min)/(n_mt-1);
                if (mt <= particle_mass[p]) continue;
                h_mt->SetBinContent(i+1, f_mt(mt, 0.0));
            }
            for (int i=0; i<n_y; ++i) {
                double y = y_min + i*(y_max-y_min)/(n_y-1);
                h_y->SetBinContent(i+1, f_y(y));
            }
            h_mt->Write();
            h_y->Write();
            delete h_mt;
            delete h_y;
        };

        if (needModel("dirac")) write_model("dirac", func_dirac_mt[p], func_dirac_y[p]);
        if (needModel("bw")) write_model("bw", func_bw_mt[p], func_bw_y[p]);
        if (needModel("ps")) write_model("ps", func_ps_mt[p], func_ps_y[p]);
    }

    f.Close();
    return 0;
}
