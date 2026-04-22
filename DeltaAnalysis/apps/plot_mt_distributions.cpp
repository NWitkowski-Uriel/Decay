// apps/plot_mt_distributions.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <getopt.h>
#include <set>
#include <sstream>

#include "TFile.h"
#include "TH1D.h"
#include "TParameter.h"
#include "TStyle.h"

#include "Constants.h"
#include "Plotting.h"

std::set<std::string> selected_dists;
std::set<std::string> selected_particles;

void parseArguments(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"distributions", required_argument, 0, 'd'},
        {"particles",     required_argument, 0, 'p'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "d:p:", long_options, nullptr)) != -1) {
        switch (c) {
            case 'd': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ','))
                    selected_dists.insert(item);
                break;
            }
            case 'p': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ','))
                    selected_particles.insert(item);
                break;
            }
        }
    }
    if (selected_dists.empty()) selected_dists = {"mt", "rapidity", "ratio"};
    if (selected_particles.empty()) {
        selected_particles = {"proton", "neutron", "piplus", "piminus", "pi0"};
    }
}

bool need(const std::string& d) { return selected_dists.count(d) > 0; }
bool needParticle(const std::string& p) { return selected_particles.count(p) > 0; }

bool hasNonZeroContent(TH1* h) {
    if (!h) return false;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        if (std::abs(h->GetBinContent(i)) > 1e-20) return true;
    }
    return false;
}

// Poprawiona funkcja ładująca dane – pomija BOM UTF-8

std::string resolveDataPath(const std::string& filename) {
    namespace fs = std::filesystem;
    if (fs::exists(filename)) return filename;
    std::string inData = "data/" + filename;
    if (fs::exists(inData)) return inData;
    return filename;
}

std::vector<std::pair<double,double>> LoadData(const std::string& filename) {
    std::vector<std::pair<double,double>> data;
    std::ifstream file(resolveDataPath(filename), std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: cannot open " << filename << std::endl;
        return data;
    }

    // Sprawdź i pomiń BOM UTF-8 (EF BB BF)
    char bom[3];
    file.read(bom, 3);
    if (!(bom[0] == (char)0xEF && bom[1] == (char)0xBB && bom[2] == (char)0xBF)) {
        file.seekg(0);
    }

    double x, y;
    int points = 0;
    while (file >> x >> y) {
        data.emplace_back(x, y);
        ++points;
    }
    file.close();
    std::cout << "Loaded " << points << " points from " << filename << std::endl;
    return data;
}

int main(int argc, char* argv[]) {
    parseArguments(argc, argv);
    if (!need("mt")) return 0;

    std::filesystem::create_directories("output");
    TFile* infile = TFile::Open("total.root", "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: cannot open total.root\n";
        return 1;
    }

    TParameter<double>* scale_param = (TParameter<double>*)infile->Get("scale_factor");
    double scale_factor = scale_param ? scale_param->GetVal() : 1.0;
    std::cout << "Scale factor = " << scale_factor << std::endl;

    struct ParticleInfo {
        std::string name; double mass; std::string xlabel; std::string title;
        std::string exp_file; int marker_style;
    };
    std::vector<ParticleInfo> all_particles = {
        {"proton",   mp,   "m_{t} - m_{p} [MeV]",   "Proton transverse mass distribution",    "protons.txt",   20},
        {"neutron",  mn,   "m_{t} - m_{n} [MeV]",   "Neutron transverse mass distribution",   "",              20},
        {"piplus",   mpi,  "m_{t} - m_{#pi} [MeV]", "#pi^{+} transverse mass distribution",   "pionpmt.txt",   21},
        {"piminus",  mpi,  "m_{t} - m_{#pi} [MeV]", "#pi^{-} transverse mass distribution",   "pionmmt.txt",   22},
        {"pi0",      mpi0, "m_{t} - m_{#pi^{0}} [MeV]", "#pi^{0} transverse mass distribution", "",              20}
    };

    for (const auto& p : all_particles) {
        if (!needParticle(p.name)) continue;

        TH1D *h_prim = (TH1D*)infile->Get(("h_total_prim_" + p.name + "_mt").c_str());
        TH1D *h_dirac = (TH1D*)infile->Get(("h_total_dirac_" + p.name + "_mt").c_str());
        TH1D *h_bw = (TH1D*)infile->Get(("h_total_bw_" + p.name + "_mt").c_str());
        TH1D *h_ps = (TH1D*)infile->Get(("h_total_ps_" + p.name + "_mt").c_str());

        if (!h_prim) {
            std::cerr << "Warning: no mt histogram for " << p.name << std::endl;
            continue;
        }

        bool has_dirac = (h_dirac && hasNonZeroContent(h_dirac));
        bool has_bw    = (h_bw && hasNonZeroContent(h_bw));
        bool has_ps    = (h_ps && hasNonZeroContent(h_ps));

        std::vector<double> x, prim, dirac, total_dirac, bw, total_bw, ps, total_ps;
        int nbins = h_prim->GetNbinsX();
        for (int i = 1; i <= nbins; ++i) {
            double xval = h_prim->GetBinCenter(i);
            if (xval - p.mass >= 0 && xval <= p.mass + 1000.0) {
                x.push_back(xval - p.mass);
                prim.push_back(h_prim->GetBinContent(i) * scale_factor);
                if (has_dirac) {
                    double d = h_dirac->GetBinContent(i) * scale_factor;
                    dirac.push_back(d);
                    total_dirac.push_back((h_prim->GetBinContent(i) + h_dirac->GetBinContent(i)) * scale_factor);
                }
                if (has_bw) {
                    double b = h_bw->GetBinContent(i) * scale_factor;
                    bw.push_back(b);
                    total_bw.push_back((h_prim->GetBinContent(i) + h_bw->GetBinContent(i)) * scale_factor);
                }
                if (has_ps) {
                    double v = h_ps->GetBinContent(i) * scale_factor;
                    ps.push_back(v);
                    total_ps.push_back((h_prim->GetBinContent(i) + h_ps->GetBinContent(i)) * scale_factor);
                }
            }
        }
        if (x.empty()) continue;

        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> series;
        series.emplace_back(&x, &prim, kRed, 1, "Primordial");
        if (has_dirac) {
            series.emplace_back(&x, &dirac, kGreen+2, 2, "From #Delta (Dirac)");
            series.emplace_back(&x, &total_dirac, kBlack, 1, "Total (Dirac)");
        }
        if (has_bw) {
            series.emplace_back(&x, &bw, kBlue, 2, "From #Delta (BW)");
            series.emplace_back(&x, &total_bw, kMagenta, 1, "Total (BW)");
        }
        if (has_ps) {
            series.emplace_back(&x, &ps, kOrange + 7, 2, "From #Delta (Phase Shift)");
            series.emplace_back(&x, &total_ps, kCyan + 2, 1, "Total (Phase Shift)");
        }

        // Dane eksperymentalne
        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> exp;
        if (!p.exp_file.empty()) {
            auto exp_data = LoadData(p.exp_file);
            if (!exp_data.empty()) {
                static std::vector<double> exp_x, exp_y;
                exp_x.clear(); exp_y.clear();
                for (const auto& pt : exp_data) {
                    if (pt.first <= 1000.0) {
                        exp_x.push_back(pt.first);
                        exp_y.push_back(pt.second);
                    }
                }
                if (!exp_x.empty()) {
                    exp.emplace_back(&exp_x, &exp_y, kBlack, p.marker_style, "Exp. data");
                }
            }
        }

        std::string filename = "output/" + p.name + "_mt_distribution";
        DrawComparisonPlot(filename, p.title, p.xlabel,
                           "1/m_{t}^{2} dN/(dm_{t} dy) [MeV^{-3}]",
                           series, exp, true, 0, 1000.0, 1e-12, 1e-4);

        if (need("ratio") && (has_dirac || has_bw || has_ps)) {
            std::vector<double> ratio_dirac, ratio_bw, ratio_ps;
            if (has_dirac) {
                for (size_t i = 0; i < x.size(); ++i) {
                    if (prim[i] > 0)
                        ratio_dirac.push_back(total_dirac[i] / prim[i]);
                    else
                        ratio_dirac.push_back(0.0);
                }
            }
            if (has_bw) {
                for (size_t i = 0; i < x.size(); ++i) {
                    if (prim[i] > 0)
                        ratio_bw.push_back(total_bw[i] / prim[i]);
                    else
                        ratio_bw.push_back(0.0);
                }
            }
            if (has_ps) {
                for (size_t i = 0; i < x.size(); ++i) {
                    if (prim[i] > 0)
                        ratio_ps.push_back(total_ps[i] / prim[i]);
                    else
                        ratio_ps.push_back(0.0);
                }
            }
            if (!ratio_dirac.empty() || !ratio_bw.empty() || !ratio_ps.empty()) {
                std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> ratio_series;
                std::vector<double> ratio_all;
                if (!ratio_dirac.empty()) {
                    ratio_series.emplace_back(&x, &ratio_dirac, kRed, 1, "Dirac");
                    ratio_all.insert(ratio_all.end(), ratio_dirac.begin(), ratio_dirac.end());
                }
                if (!ratio_bw.empty()) {
                    ratio_series.emplace_back(&x, &ratio_bw, kBlue, 1, "Breit‑Wigner");
                    ratio_all.insert(ratio_all.end(), ratio_bw.begin(), ratio_bw.end());
                }
                if (!ratio_ps.empty()) {
                    ratio_series.emplace_back(&x, &ratio_ps, kOrange + 7, 1, "Phase Shift");
                    ratio_all.insert(ratio_all.end(), ratio_ps.begin(), ratio_ps.end());
                }
                double ymin = 0.9 * *std::min_element(ratio_all.begin(), ratio_all.end());
                double ymax = 1.1 * *std::max_element(ratio_all.begin(), ratio_all.end());
                DrawComparisonPlot(filename + "_ratio",
                                   p.title + " (total/primordial ratio)",
                                   p.xlabel, "total / primordial",
                                   ratio_series, {}, false, 0, 1000.0, ymin, ymax);
            }
        }
    }

    infile->Close();
    return 0;
}