// compute_delta_decays.cpp
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include "TFile.h"
#include "TH1D.h"
#include "Constants.h"
#include "DecayFunctions.h"

struct Args {
    std::string particles; // lista oddzielona przecinkami
    std::string model;     // "dirac", "bw", "both"
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
        }
    }
    if (args.particles.empty()) args.particles = "proton,neutron,piplus,piminus,pi0";
    if (args.model.empty()) args.model = "both";
}

bool needParticle(const std::string& name, const std::string& list) {
    return list.find(name) != std::string::npos;
}

int main(int argc, char* argv[]) {
    parse(argc, argv);
    const int n_mt = 1000;
    double mt_min = 0.0, mt_max = 3000.0;
    const int n_y = 101;
    double y_min = -2.0, y_max = 2.0;

    TFile f("delta.root", "RECREATE");

    // Funkcje do wywołania dla danej cząstki i modelu
    std::map<std::string, std::function<double(double,double)>> func_dirac_mt;
    std::map<std::string, std::function<double(double,double)>> func_bw_mt;
    std::map<std::string, std::function<double(double)>> func_dirac_y;
    std::map<std::string, std::function<double(double)>> func_bw_y;

    func_dirac_mt["proton"]   = dN_dmt_proton_from_Delta_Dirac;
    func_bw_mt["proton"]      = dN_dmt_proton_from_Delta_BW;
    func_dirac_y["proton"]    = [](double y) { return dN_dy_proton_from_Delta_Dirac(y); }; // trzeba zaimplementować w DecayFunctions
    func_bw_y["proton"]       = [](double y) { return dN_dy_proton_from_Delta_BW(y); };
    // analogicznie dla neutron, piplus, piminus, pi0...

    std::vector<std::string> particles = {"proton","neutron","piplus","piminus","pi0"};
    for (const auto& p : particles) {
        if (!needParticle(p, args.particles)) continue;

        if (args.model == "dirac" || args.model == "both") {
            TH1D *h_mt = new TH1D(("h_delta_dirac_" + p + "_mt").c_str(), "", n_mt, mt_min, mt_max);
            TH1D *h_y  = new TH1D(("h_delta_dirac_" + p + "_y").c_str(),  "", n_y, y_min, y_max);
            for (int i=0; i<n_mt; ++i) {
                double mt = mt_min + i*(mt_max-mt_min)/(n_mt-1);
                if (mt <= mp && p=="proton") continue;
                h_mt->SetBinContent(i+1, func_dirac_mt[p](mt, 0.0));
            }
            for (int i=0; i<n_y; ++i) {
                double y = y_min + i*(y_max-y_min)/(n_y-1);
                h_y->SetBinContent(i+1, func_dirac_y[p](y));
            }
            h_mt->Write(); h_y->Write();
            delete h_mt; delete h_y;
        }
        if (args.model == "bw" || args.model == "both") {
            TH1D *h_mt = new TH1D(("h_delta_bw_" + p + "_mt").c_str(), "", n_mt, mt_min, mt_max);
            TH1D *h_y  = new TH1D(("h_delta_bw_" + p + "_y").c_str(),  "", n_y, y_min, y_max);
            // podobnie wypełnianie
            h_mt->Write(); h_y->Write();
        }
    }
    f.Close();
    return 0;
}