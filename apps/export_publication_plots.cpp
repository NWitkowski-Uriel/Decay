#include <algorithm>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "TColor.h"
#include "Constants.h"
#include "DecayFunctions.h"
#include "Plotting.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

struct Args {
    std::set<std::string> particles;
    std::string output_prefix = "output/publication";
    int n_mt = 500;
    double mt_max_offset = 1000.0;
    bool include_primordial = true;
} args;

void parse_args(int argc, char* argv[]) {
    static option long_options[] = {
        {"particles", required_argument, nullptr, 'p'},
        {"output-prefix", required_argument, nullptr, 'o'},
        {"nmt", required_argument, nullptr, 'n'},
        {"mt-max-offset", required_argument, nullptr, 'm'},
        {"no-primordial", no_argument, nullptr, 'x'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "p:o:n:m:x", long_options, nullptr)) != -1) {
        switch (c) {
            case 'p': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    if (!item.empty()) args.particles.insert(item);
                }
                break;
            }
            case 'o': args.output_prefix = optarg; break;
            case 'n': args.n_mt = std::max(100, std::atoi(optarg)); break;
            case 'm': args.mt_max_offset = std::max(100.0, std::atof(optarg)); break;
            case 'x': args.include_primordial = false; break;
            default: break;
        }
    }

    if (args.particles.empty()) {
        args.particles = {"proton", "neutron", "piplus", "piminus", "pi0"};
    }
}

struct ParticleInfo {
    double mass;
    double mu;
    double g;
    std::string title;
    std::string x_title;
    std::function<double(double,double)> dirac_mt;
    std::function<double(double,double)> bw_mt;
    std::function<double(double,double)> ps_mt;
};

std::map<std::string, ParticleInfo> build_particle_map() {
    std::map<std::string, ParticleInfo> m;
    m["proton"] = {mp, mu_p, g_proton,
        "Proton transverse-mass spectrum", "m_{T} - m_{p} [MeV]",
        dN_dmt_proton_from_Delta_Dirac, dN_dmt_proton_from_Delta_BW, dN_dmt_proton_from_Delta_PS};
    m["neutron"] = {mn, mu_n, g_neutron,
        "Neutron transverse-mass spectrum", "m_{T} - m_{n} [MeV]",
        dN_dmt_neutron_from_Delta_Dirac, dN_dmt_neutron_from_Delta_BW, dN_dmt_neutron_from_Delta_PS};
    m["piplus"] = {mpi, mu_piplus, g_pion,
        "#pi^{+} transverse-mass spectrum", "m_{T} - m_{#pi^{+}} [MeV]",
        dN_dmt_piplus_from_Delta_Dirac, dN_dmt_piplus_from_Delta_BW, dN_dmt_piplus_from_Delta_PS};
    m["piminus"] = {mpi, mu_piminus, g_pion,
        "#pi^{-} transverse-mass spectrum", "m_{T} - m_{#pi^{-}} [MeV]",
        dN_dmt_piminus_from_Delta_Dirac, dN_dmt_piminus_from_Delta_BW, dN_dmt_piminus_from_Delta_PS};
    m["pi0"] = {mpi0, mu_pi0, g_pion,
        "#pi^{0} transverse-mass spectrum", "m_{T} - m_{#pi^{0}} [MeV]",
        dN_dmt_pi0_from_Delta_Dirac, dN_dmt_pi0_from_Delta_BW, dN_dmt_pi0_from_Delta_PS};
    return m;
}

void build_x_grid(double mass, std::vector<double>& x_raw, std::vector<double>& x_shifted) {
    x_raw.clear();
    x_shifted.clear();
    const double mt_min = mass + 1e-3;
    const double mt_max = mass + args.mt_max_offset;
    const double step = (mt_max - mt_min) / (args.n_mt - 1);
    for (int i = 0; i < args.n_mt; ++i) {
        const double mt = mt_min + i * step;
        x_raw.push_back(mt);
        x_shifted.push_back(mt - mass);
    }
}

std::vector<double> make_ratio(const std::vector<double>& num, const std::vector<double>& den) {
    std::vector<double> r(num.size(), 0.0);
    for (size_t i = 0; i < num.size(); ++i) {
        r[i] = (i < den.size() && den[i] > 0.0) ? num[i] / den[i] : 0.0;
    }
    return r;
}

} // namespace

int main(int argc, char* argv[]) {
    parse_args(argc, argv);
    RadialGrid::instance().setNumPoints(200);
    const auto particles = build_particle_map();

    for (const auto& name : args.particles) {
        auto it = particles.find(name);
        if (it == particles.end()) continue;
        const auto& p = it->second;

        std::vector<double> mt_raw, mt_plot;
        build_x_grid(p.mass, mt_raw, mt_plot);

        std::vector<double> prim(args.n_mt), dirac(args.n_mt), bw(args.n_mt), ps(args.n_mt);
        std::vector<double> total_dirac(args.n_mt), total_bw(args.n_mt), total_ps(args.n_mt);

        for (int i = 0; i < args.n_mt; ++i) {
            const double mt = mt_raw[i];
            prim[i] = dN_dmt_primordial(mt, 0.0, p.mu, p.mass, p.g);
            dirac[i] = p.dirac_mt(mt, 0.0);
            bw[i] = p.bw_mt(mt, 0.0);
            ps[i] = p.ps_mt(mt, 0.0);
            total_dirac[i] = prim[i] + dirac[i];
            total_bw[i] = prim[i] + bw[i];
            total_ps[i] = prim[i] + ps[i];
        }

        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> spectra_series;
        if (args.include_primordial) spectra_series.emplace_back(&mt_plot, &prim, kBlack, 1, "Primordial");
        spectra_series.emplace_back(&mt_plot, &total_dirac, kRed + 1, 1, "Total Dirac");
        spectra_series.emplace_back(&mt_plot, &total_bw, kBlue + 1, 1, "Total BW");
        spectra_series.emplace_back(&mt_plot, &total_ps, kGreen + 2, 1, "Total PS");

        DrawComparisonPlot(
            args.output_prefix + "_" + name + "_spectra",
            p.title,
            p.x_title,
            "1/m_{T}^{2} dN/(dm_{T} dy) [MeV^{-3}]",
            spectra_series,
            {},
            true,
            0.0,
            args.mt_max_offset,
            1e-12,
            1e-2
        );

        std::vector<double> ratio_bw_dirac = make_ratio(total_bw, total_dirac);
        std::vector<double> ratio_ps_dirac = make_ratio(total_ps, total_dirac);
        std::vector<double> ratio_ps_bw = make_ratio(total_ps, total_bw);

        std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>> ratio_series;
        ratio_series.emplace_back(&mt_plot, &ratio_bw_dirac, kBlue + 1, 1, "BW / Dirac");
        ratio_series.emplace_back(&mt_plot, &ratio_ps_dirac, kGreen + 2, 1, "PS / Dirac");
        ratio_series.emplace_back(&mt_plot, &ratio_ps_bw, kMagenta + 1, 2, "PS / BW");

        DrawComparisonPlot(
            args.output_prefix + "_" + name + "_ratios",
            p.title + " ratios",
            p.x_title,
            "Model ratio",
            ratio_series,
            {},
            false,
            0.0,
            args.mt_max_offset,
            0.5,
            1.5
        );
    }

    std::cout << "Exported publication plots with prefix: " << args.output_prefix << std::endl;
    return 0;
}
