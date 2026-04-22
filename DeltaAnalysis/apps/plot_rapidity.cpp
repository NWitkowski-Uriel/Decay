// apps/plot_rapidity.cpp (uproszczona, bez odbicia)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
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


std::string resolveDataPath(const std::string& filename) {
    namespace fs = std::filesystem;
    if (fs::exists(filename)) return filename;
    std::string inData = "data/" + filename;
    if (fs::exists(inData)) return inData;
    return filename;
}

std::vector<std::pair<double,double>> LoadData(const std::string& filename) {
    std::vector<std::pair<double,double>> data;
    std::ifstream file(resolveDataPath(filename));
    if (!file.is_open()) {
        std::cerr << "Warning: cannot open " << filename << std::endl;
        return data;
    }
    double x, y;
    while (file >> x >> y) data.emplace_back(x, y);
    return data;
}

int main(int argc, char* argv[]) {
    parseArguments(argc, argv);
    if (!need("rapidity") && !need("ratio")) return 0;

    std::filesystem::create_directories("output");
    TFile* infile = TFile::Open("total.root", "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: cannot open total.root\n";
        return 1;
    }

    TParameter<double>* scale_param = (TParameter<double>*)infile->Get("scale_factor");
    double scale_factor = scale_param ? scale_param->GetVal() : 1.0;

    struct ParticleInfo {
        std::string name; std::string title; std::string exp_file; int marker_style;
    };
    std::vector<ParticleInfo> all_particles = {
        {"proton",   "Proton rapidity distribution",    "protonsY.txt",   20},
        {"neutron",  "Neutron rapidity distribution",   "",               20},
        {"piplus",   "#pi^{+} rapidity distribution",   "pionpY.txt",     21},
        {"piminus",  "#pi^{-} rapidity distribution",   "pionmY.txt",     22},
        {"pi0",      "#pi^{0} rapidity distribution",   "",               20}
    };

    for (const auto& p : all_particles) {
        if (!needParticle(p.name)) continue;

        TH1D *h_prim = (TH1D*)infile->Get(("h_total_prim_" + p.name + "_y").c_str());
        TH1D *h_dirac = (TH1D*)infile->Get(("h_total_dirac_" + p.name + "_y").c_str());
        TH1D *h_bw = (TH1D*)infile->Get(("h_total_bw_" + p.name + "_y").c_str());
        TH1D *h_ps = (TH1D*)infile->Get(("h_total_ps_" + p.name + "_y").c_str());

        if (!h_prim) {
            std::cerr << "Warning: no y histogram for " << p.name << std::endl;
            continue;
        }

        bool has_dirac = (h_dirac && hasNonZeroContent(h_dirac));
        bool has_bw    = (h_bw && hasNonZeroContent(h_bw));
        bool has_ps    = (h_ps && hasNonZeroContent(h_ps));

        // Odczytaj pełny zakres y (bez odbicia)
        int nbins = h_prim->GetNbinsX();
        std::vector<double> y_full, prim_full, dirac_full, bw_full, ps_full, total_dirac_full, total_bw_full, total_ps_full;
        for (int i = 1; i <= nbins; ++i) {
            double y = h_prim->GetBinCenter(i);
            y_full.push_back(y);
            prim_full.push_back(h_prim->GetBinContent(i) * scale_factor);
            if (has_dirac) {
                double d = h_dirac->GetBinContent(i) * scale_factor;
                dirac_full.push_back(d);
                total_dirac_full.push_back((h_prim->GetBinContent(i) + h_dirac->GetBinContent(i)) * scale_factor);
            }
            if (has_bw) {
                double b = h_bw->GetBinContent(i) * scale_factor;
                bw_full.push_back(b);
                total_bw_full.push_back((h_prim->GetBinContent(i) + h_bw->GetBinContent(i)) * scale_factor);
            }
            if (has_ps) {
                double ps = h_ps->GetBinContent(i) * scale_factor;
                ps_full.push_back(ps);
                total_ps_full.push_back((h_prim->GetBinContent(i) + h_ps->GetBinContent(i)) * scale_factor);
            }
        }

        // Serie modelowe
        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> series;
        series.emplace_back(&y_full, &prim_full, kRed, 1, "Primordial");
        if (has_dirac) {
            series.emplace_back(&y_full, &dirac_full, kGreen+2, 2, "From #Delta (Dirac)");
            series.emplace_back(&y_full, &total_dirac_full, kBlack, 1, "Total (Dirac)");
        }
        if (has_bw) {
            series.emplace_back(&y_full, &bw_full, kBlue, 2, "From #Delta (BW)");
            series.emplace_back(&y_full, &total_bw_full, kMagenta, 1, "Total (BW)");
        }
        if (has_ps) {
            series.emplace_back(&y_full, &ps_full, kOrange + 7, 2, "From #Delta (Phase Shift)");
            series.emplace_back(&y_full, &total_ps_full, kCyan + 2, 1, "Total (Phase Shift)");
        }

        // Dane eksperymentalne
        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> exp;
        if (!p.exp_file.empty()) {
            auto exp_data = LoadData(p.exp_file);
            if (!exp_data.empty()) {
                static std::vector<double> exp_x, exp_y;
                exp_x.clear(); exp_y.clear();
                for (const auto& pt : exp_data) {
                    exp_x.push_back(pt.first);
                    exp_y.push_back(pt.second);
                }
                if (!exp_x.empty()) {
                    exp.emplace_back(&exp_x, &exp_y, kBlack, p.marker_style, "Exp. data");
                }
            }
        }

        if (need("rapidity")) {
            double max_y = 0.0;
            if (has_dirac) max_y = *std::max_element(total_dirac_full.begin(), total_dirac_full.end());
            else if (has_bw) max_y = *std::max_element(total_bw_full.begin(), total_bw_full.end());
            else if (has_ps) max_y = *std::max_element(total_ps_full.begin(), total_ps_full.end());
            else max_y = *std::max_element(prim_full.begin(), prim_full.end());
            if (max_y <= 0) max_y = 1.0;
            DrawComparisonPlot("output/" + p.name + "_rapidity",
                               p.title, "y", "dN/dy",
                               series, exp, false, -2.0, 2.0, 0, max_y * 1.2);
        }

        if (need("ratio") && (has_dirac || has_bw || has_ps)) {
            std::vector<double> ratio_dirac(y_full.size()), ratio_bw(y_full.size()), ratio_ps(y_full.size());
            if (has_dirac) {
                for (size_t i = 0; i < y_full.size(); ++i) {
                    if (prim_full[i] > 0)
                        ratio_dirac[i] = total_dirac_full[i] / prim_full[i];
                    else
                        ratio_dirac[i] = 0.0;
                }
            }
            if (has_bw) {
                for (size_t i = 0; i < y_full.size(); ++i) {
                    if (prim_full[i] > 0)
                        ratio_bw[i] = total_bw_full[i] / prim_full[i];
                    else
                        ratio_bw[i] = 0.0;
                }
            }
            if (has_ps) {
                for (size_t i = 0; i < y_full.size(); ++i) {
                    if (prim_full[i] > 0)
                        ratio_ps[i] = total_ps_full[i] / prim_full[i];
                    else
                        ratio_ps[i] = 0.0;
                }
            }
            std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> ratio_series;
            std::vector<double> ratio_all;
            if (has_dirac) {
                ratio_series.emplace_back(&y_full, &ratio_dirac, kRed, 1, "Dirac");
                ratio_all.insert(ratio_all.end(), ratio_dirac.begin(), ratio_dirac.end());
            }
            if (has_bw) {
                ratio_series.emplace_back(&y_full, &ratio_bw, kBlue, 1, "Breit‑Wigner");
                ratio_all.insert(ratio_all.end(), ratio_bw.begin(), ratio_bw.end());
            }
            if (has_ps) {
                ratio_series.emplace_back(&y_full, &ratio_ps, kOrange + 7, 1, "Phase Shift");
                ratio_all.insert(ratio_all.end(), ratio_ps.begin(), ratio_ps.end());
            }
            if (!ratio_series.empty()) {
                double ymin = 0.9 * *std::min_element(ratio_all.begin(), ratio_all.end());
                double ymax = 1.1 * *std::max_element(ratio_all.begin(), ratio_all.end());
                DrawComparisonPlot("output/" + p.name + "_rapidity_ratio",
                                   p.title + " (total/primordial ratio)",
                                   "y", "total / primordial",
                                   ratio_series, {}, false, -2.0, 2.0, ymin, ymax);
            }
        }
    }

    infile->Close();
    return 0;
}