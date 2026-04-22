// apps/compute_yields.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <getopt.h>
#include <set>
#include <sstream>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TStopwatch.h"

#include "Constants.h"
#include "MathUtils.h"
#include "Thermodynamics.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"

// ============================================================================
// Struktura do przechowywania wyboru użytkownika
// ============================================================================
struct Selection {
    std::set<std::string> particles;   // "proton","neutron","piplus","piminus","pi0"
    std::set<std::string> dists;       // "mt","rapidity","ratio"
    std::string model;                 // "dirac","bw","both"
    std::string output;                // nazwa pliku wynikowego
} sel;

// ============================================================================
// Parsowanie argumentów linii poleceń
// ============================================================================
void parseArguments(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"particles",     required_argument, 0, 'p'},
        {"distributions", required_argument, 0, 'd'},
        {"model",         required_argument, 0, 'm'},
        {"output",        required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "p:d:m:o:", long_options, nullptr)) != -1) {
        switch (c) {
            case 'p': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    sel.particles.insert(item);
                }
                break;
            }
            case 'd': {
                std::stringstream ss(optarg);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    sel.dists.insert(item);
                }
                break;
            }
            case 'm': sel.model = optarg; break;
            case 'o': sel.output = optarg; break;
            default: break;
        }
    }
    // Domyślne wartości
    if (sel.particles.empty()) {
        sel.particles = {"proton","neutron","piplus","piminus","pi0"};
    }
    if (sel.dists.empty()) {
        sel.dists = {"mt","rapidity","ratio"};
    }
    if (sel.model.empty()) sel.model = "both";
    if (sel.output.empty()) sel.output = "output/yields.root";
}

// ============================================================================
// Funkcje pomocnicze
// ============================================================================
bool need(const std::string& dist) {
    return sel.dists.count(dist) > 0;
}

bool needParticle(const std::string& p) {
    return sel.particles.count(p) > 0;
}

bool needModel(const std::string& m) {
    if (sel.model == "both") return true;
    return sel.model == m;
}

// ============================================================================
// MAIN
// ============================================================================
int main(int argc, char* argv[]) {
    parseArguments(argc, argv);

    // Utwórz katalog dla pliku wyjściowego, jeśli nie istnieje
    std::filesystem::path out_path(sel.output);
    std::filesystem::create_directories(out_path.parent_path());

    // Inicjalizacja siatki radialnej (200 punktów)
    RadialGrid::instance().setNumPoints(200);

    // Konfiguracja OpenMP
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    std::cout << "==========================================" << std::endl;
    std::cout << "COMPUTE YIELDS – SELECTIVE CALCULATION" << std::endl;
    std::cout << "Using " << max_threads << " threads" << std::endl;
    std::cout << "Selected particles: ";
    for (const auto& p : sel.particles) std::cout << p << " ";
    std::cout << std::endl;
    std::cout << "Selected distributions: ";
    for (const auto& d : sel.dists) std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "Model: " << sel.model << std::endl;
    std::cout << "Output file: " << sel.output << std::endl;
    std::cout << "==========================================" << std::endl;

    TStopwatch timer;
    timer.Start();

    // Siatki dla mt i y (używamy symetrii dla y)
    const int npoints_mt = 1000;
    const double mt_min = 0.0;
    const double mt_max = 3000.0;
    const double dmt = (mt_max - mt_min) / (npoints_mt - 1);

    const int npoints_y_half = 101;
    const double y_min_half = 0.0;
    const double y_max_half = 2.0;
    const double dy = (y_max_half - y_min_half) / (npoints_y_half - 1);
    const int npoints_y_full = 2 * npoints_y_half - 1;

    std::vector<double> mt_vals(npoints_mt);
    for (int i = 0; i < npoints_mt; ++i) mt_vals[i] = mt_min + i * dmt;

    std::vector<double> y_half(npoints_y_half);
    for (int i = 0; i < npoints_y_half; ++i) y_half[i] = y_min_half + i * dy;

    // Struktury do przechowywania wyników (wektory)
    struct MtData {
        std::vector<double> prim, dirac, bw, total_dirac, total_bw;
    };
    struct YData {
        std::vector<double> prim, dirac, bw, total_dirac, total_bw;
    };
    std::map<std::string, MtData> mt_data;
    std::map<std::string, YData> y_data;

    // Wektory dla delt (potrzebne do rapidity)
    std::vector<double> dndy_Dpp_Dirac, dndy_Dp_Dirac, dndy_D0_Dirac, dndy_Dm_Dirac;
    std::vector<double> dndy_Dpp_BW, dndy_Dp_BW, dndy_D0_BW, dndy_Dm_BW;
    if (need("rapidity")) {
        dndy_Dpp_Dirac.resize(npoints_y_half);
        dndy_Dp_Dirac.resize(npoints_y_half);
        dndy_D0_Dirac.resize(npoints_y_half);
        dndy_Dm_Dirac.resize(npoints_y_half);
        dndy_Dpp_BW.resize(npoints_y_half);
        dndy_Dp_BW.resize(npoints_y_half);
        dndy_D0_BW.resize(npoints_y_half);
        dndy_Dm_BW.resize(npoints_y_half);
    }

    // Inicjalizacja wektorów dla wybranych cząstek
    for (const auto& p : sel.particles) {
        if (need("mt")) {
            mt_data[p].prim.resize(npoints_mt, 0.0);
            mt_data[p].dirac.resize(npoints_mt, 0.0);
            mt_data[p].bw.resize(npoints_mt, 0.0);
            mt_data[p].total_dirac.resize(npoints_mt, 0.0);
            mt_data[p].total_bw.resize(npoints_mt, 0.0);
        }
        if (need("rapidity")) {
            y_data[p].prim.resize(npoints_y_half, 0.0);
            y_data[p].dirac.resize(npoints_y_half, 0.0);
            y_data[p].bw.resize(npoints_y_half, 0.0);
            y_data[p].total_dirac.resize(npoints_y_half, 0.0);
            y_data[p].total_bw.resize(npoints_y_half, 0.0);
        }
    }

    // -------------------------------------------------------------------------
    // Obliczenia widm poprzecznych (mt) dla y=0
    // -------------------------------------------------------------------------
    if (need("mt")) {
        std::cout << "\nComputing transverse mass spectra (y=0) ..." << std::endl;
        #pragma omp parallel for
        for (int i = 0; i < npoints_mt; ++i) {
            double mt = mt_vals[i];
            if (i % 200 == 0) {
                #pragma omp critical
                std::cout << "  mt index " << i << " / " << npoints_mt << std::endl;
            }

            if (needParticle("proton") && mt > mp) {
                double prim = dN_dmt_primordial(mt, 0.0, mu_p, mp, g_proton);
                mt_data["proton"].prim[i] = prim;
                if (needModel("dirac")) {
                    double dirac = dN_dmt_proton_from_Delta_Dirac(mt, 0.0);
                    mt_data["proton"].dirac[i] = dirac;
                    mt_data["proton"].total_dirac[i] = prim + dirac;
                }
                if (needModel("bw")) {
                    double bw = dN_dmt_proton_from_Delta_BW(mt, 0.0);
                    mt_data["proton"].bw[i] = bw;
                    mt_data["proton"].total_bw[i] = prim + bw;
                }
            }
            if (needParticle("neutron") && mt > mn) {
                double prim = dN_dmt_primordial(mt, 0.0, mu_n, mn, g_neutron);
                mt_data["neutron"].prim[i] = prim;
                if (needModel("dirac")) {
                    double dirac = dN_dmt_neutron_from_Delta_Dirac(mt, 0.0);
                    mt_data["neutron"].dirac[i] = dirac;
                    mt_data["neutron"].total_dirac[i] = prim + dirac;
                }
                if (needModel("bw")) {
                    double bw = dN_dmt_neutron_from_Delta_BW(mt, 0.0);
                    mt_data["neutron"].bw[i] = bw;
                    mt_data["neutron"].total_bw[i] = prim + bw;
                }
            }
            if (needParticle("piplus") && mt > mpi) {
                double prim = dN_dmt_primordial(mt, 0.0, mu_piplus, mpi, g_pion);
                mt_data["piplus"].prim[i] = prim;
                if (needModel("dirac")) {
                    double dirac = dN_dmt_piplus_from_Delta_Dirac(mt, 0.0);
                    mt_data["piplus"].dirac[i] = dirac;
                    mt_data["piplus"].total_dirac[i] = prim + dirac;
                }
                if (needModel("bw")) {
                    double bw = dN_dmt_piplus_from_Delta_BW(mt, 0.0);
                    mt_data["piplus"].bw[i] = bw;
                    mt_data["piplus"].total_bw[i] = prim + bw;
                }
            }
            if (needParticle("piminus") && mt > mpi) {
                double prim = dN_dmt_primordial(mt, 0.0, mu_piminus, mpi, g_pion);
                mt_data["piminus"].prim[i] = prim;
                if (needModel("dirac")) {
                    double dirac = dN_dmt_piminus_from_Delta_Dirac(mt, 0.0);
                    mt_data["piminus"].dirac[i] = dirac;
                    mt_data["piminus"].total_dirac[i] = prim + dirac;
                }
                if (needModel("bw")) {
                    double bw = dN_dmt_piminus_from_Delta_BW(mt, 0.0);
                    mt_data["piminus"].bw[i] = bw;
                    mt_data["piminus"].total_bw[i] = prim + bw;
                }
            }
            if (needParticle("pi0") && mt > mpi0) {
                double prim = dN_dmt_primordial(mt, 0.0, mu_pi0, mpi0, g_pion);
                mt_data["pi0"].prim[i] = prim;
                if (needModel("dirac")) {
                    double dirac = dN_dmt_pi0_from_Delta_Dirac(mt, 0.0);
                    mt_data["pi0"].dirac[i] = dirac;
                    mt_data["pi0"].total_dirac[i] = prim + dirac;
                }
                if (needModel("bw")) {
                    double bw = dN_dmt_pi0_from_Delta_BW(mt, 0.0);
                    mt_data["pi0"].bw[i] = bw;
                    mt_data["pi0"].total_bw[i] = prim + bw;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Obliczenia rozkładów rapidity (y>=0)
    // -------------------------------------------------------------------------
    if (need("rapidity")) {
        std::cout << "\nComputing rapidity distributions (y>=0) ..." << std::endl;
        #pragma omp parallel for
        for (int i = 0; i < npoints_y_half; ++i) {
            double y = y_half[i];
            if (i % 20 == 0) {
                #pragma omp critical
                std::cout << "  y index " << i << " / " << npoints_y_half << " (y = " << y << ")" << std::endl;
            }

            if (needModel("dirac")) {
                dndy_Dpp_Dirac[i] = dN_dy_Delta_Dirac(y, mu_Dpp, g_Delta);
                dndy_Dp_Dirac[i]  = dN_dy_Delta_Dirac(y, mu_Dp,  g_Delta);
                dndy_D0_Dirac[i]  = dN_dy_Delta_Dirac(y, mu_D0,  g_Delta);
                dndy_Dm_Dirac[i]  = dN_dy_Delta_Dirac(y, mu_Dm,  g_Delta);
            }
            if (needModel("bw")) {
                dndy_Dpp_BW[i] = dN_dy_Delta_BW(y, mu_Dpp, g_Delta);
                dndy_Dp_BW[i]  = dN_dy_Delta_BW(y, mu_Dp,  g_Delta);
                dndy_D0_BW[i]  = dN_dy_Delta_BW(y, mu_D0,  g_Delta);
                dndy_Dm_BW[i]  = dN_dy_Delta_BW(y, mu_Dm,  g_Delta);
            }

            for (const auto& p : sel.particles) {
                double prim = 0.0;
                if (p == "proton")   prim = dN_dy_primordial_full(y, mu_p, mp, g_proton);
                if (p == "neutron")  prim = dN_dy_primordial_full(y, mu_n, mn, g_neutron);
                if (p == "piplus")   prim = dN_dy_primordial_full(y, mu_piplus, mpi, g_pion);
                if (p == "piminus")  prim = dN_dy_primordial_full(y, mu_piminus, mpi, g_pion);
                if (p == "pi0")      prim = dN_dy_primordial_full(y, mu_pi0, mpi0, g_pion);
                y_data[p].prim[i] = prim;

                if (needModel("dirac")) {
                    double dirac = 0.0;
                    if (p == "proton")   dirac = dndy_Dpp_Dirac[i] + (2.0/3.0)*dndy_Dp_Dirac[i] + (1.0/3.0)*dndy_D0_Dirac[i];
                    if (p == "neutron")  dirac = (1.0/3.0)*dndy_Dp_Dirac[i] + (2.0/3.0)*dndy_D0_Dirac[i] + dndy_Dm_Dirac[i];
                    if (p == "piplus")   dirac = dndy_Dpp_Dirac[i] + (1.0/3.0)*dndy_Dp_Dirac[i];
                    if (p == "piminus")  dirac = (1.0/3.0)*dndy_D0_Dirac[i] + dndy_Dm_Dirac[i];
                    if (p == "pi0")      dirac = (2.0/3.0)*dndy_Dp_Dirac[i] + (2.0/3.0)*dndy_D0_Dirac[i];
                    y_data[p].dirac[i] = dirac;
                    y_data[p].total_dirac[i] = prim + dirac;
                }
                if (needModel("bw")) {
                    double bw = 0.0;
                    if (p == "proton")   bw = dndy_Dpp_BW[i] + (2.0/3.0)*dndy_Dp_BW[i] + (1.0/3.0)*dndy_D0_BW[i];
                    if (p == "neutron")  bw = (1.0/3.0)*dndy_Dp_BW[i] + (2.0/3.0)*dndy_D0_BW[i] + dndy_Dm_BW[i];
                    if (p == "piplus")   bw = dndy_Dpp_BW[i] + (1.0/3.0)*dndy_Dp_BW[i];
                    if (p == "piminus")  bw = (1.0/3.0)*dndy_D0_BW[i] + dndy_Dm_BW[i];
                    if (p == "pi0")      bw = (2.0/3.0)*dndy_Dp_BW[i] + (2.0/3.0)*dndy_D0_BW[i];
                    y_data[p].bw[i] = bw;
                    y_data[p].total_bw[i] = prim + bw;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Obliczenia całkowitych yieldów protonu (dla skalowania)
    // -------------------------------------------------------------------------
    double total_proton_prim = 0.0, total_proton_dirac = 0.0, total_proton_bw = 0.0;
    if (need("rapidity") && needParticle("proton")) {
        double dy_full = (y_max_half * 2) / (npoints_y_full - 1);
        for (int i = 0; i < npoints_y_half; ++i) {
            double weight = (i == 0) ? 0.5 : 1.0;
            total_proton_prim += weight * y_data["proton"].prim[i] * dy_full;
            if (needModel("dirac")) total_proton_dirac += weight * y_data["proton"].dirac[i] * dy_full;
            if (needModel("bw"))    total_proton_bw += weight * y_data["proton"].bw[i] * dy_full;
        }
        total_proton_prim = 2 * total_proton_prim - (y_data["proton"].prim[0] * dy_full) / 2;
        if (needModel("dirac")) total_proton_dirac = 2 * total_proton_dirac - (y_data["proton"].dirac[0] * dy_full) / 2;
        if (needModel("bw"))    total_proton_bw = 2 * total_proton_bw - (y_data["proton"].bw[0] * dy_full) / 2;
    }

    double total_proton_total_dirac = total_proton_prim + total_proton_dirac;
    double total_proton_total_bw    = total_proton_prim + total_proton_bw;
    double scale_factor = (total_proton_total_dirac > 0) ? N_proton_exp / total_proton_total_dirac : 1.0;

    // -------------------------------------------------------------------------
    // Zapis do pliku ROOT – poprawiony sposób
    // -------------------------------------------------------------------------
    std::cout << "\nSaving results to " << sel.output << " ..." << std::endl;
    TFile* outfile = TFile::Open(sel.output.c_str(), "RECREATE");

    // -------------------------------------------------------------------------
    // Drzewo dla widm poprzecznych (mt)
    // -------------------------------------------------------------------------
    TTree* mt_tree = nullptr;
    if (need("mt")) {
        mt_tree = new TTree("mt_tree", "Transverse mass spectra at y=0");

        // Zmienne do przechowywania wartości dla każdej cząstki
        double mt, mt_minus_mp, mt_minus_mn, mt_minus_mpi, mt_minus_mpi0;
        mt_tree->Branch("mt", &mt, "mt/D");
        mt_tree->Branch("mt_minus_mp", &mt_minus_mp, "mt_minus_mp/D");
        mt_tree->Branch("mt_minus_mn", &mt_minus_mn, "mt_minus_mn/D");
        mt_tree->Branch("mt_minus_mpi", &mt_minus_mpi, "mt_minus_mpi/D");
        mt_tree->Branch("mt_minus_mpi0", &mt_minus_mpi0, "mt_minus_mpi0/D");

        // Dla każdej cząstki tworzymy gałęzie
        if (needParticle("proton")) {
            mt_tree->Branch("proton_prim", &mt_data["proton"].prim[0]);
            if (needModel("dirac")) {
                mt_tree->Branch("proton_dirac", &mt_data["proton"].dirac[0]);
                mt_tree->Branch("proton_total_dirac", &mt_data["proton"].total_dirac[0]);
            }
            if (needModel("bw")) {
                mt_tree->Branch("proton_bw", &mt_data["proton"].bw[0]);
                mt_tree->Branch("proton_total_bw", &mt_data["proton"].total_bw[0]);
            }
        }
        if (needParticle("neutron")) {
            mt_tree->Branch("neutron_prim", &mt_data["neutron"].prim[0]);
            if (needModel("dirac")) {
                mt_tree->Branch("neutron_dirac", &mt_data["neutron"].dirac[0]);
                mt_tree->Branch("neutron_total_dirac", &mt_data["neutron"].total_dirac[0]);
            }
            if (needModel("bw")) {
                mt_tree->Branch("neutron_bw", &mt_data["neutron"].bw[0]);
                mt_tree->Branch("neutron_total_bw", &mt_data["neutron"].total_bw[0]);
            }
        }
        if (needParticle("piplus")) {
            mt_tree->Branch("piplus_prim", &mt_data["piplus"].prim[0]);
            if (needModel("dirac")) {
                mt_tree->Branch("piplus_dirac", &mt_data["piplus"].dirac[0]);
                mt_tree->Branch("piplus_total_dirac", &mt_data["piplus"].total_dirac[0]);
            }
            if (needModel("bw")) {
                mt_tree->Branch("piplus_bw", &mt_data["piplus"].bw[0]);
                mt_tree->Branch("piplus_total_bw", &mt_data["piplus"].total_bw[0]);
            }
        }
        if (needParticle("piminus")) {
            mt_tree->Branch("piminus_prim", &mt_data["piminus"].prim[0]);
            if (needModel("dirac")) {
                mt_tree->Branch("piminus_dirac", &mt_data["piminus"].dirac[0]);
                mt_tree->Branch("piminus_total_dirac", &mt_data["piminus"].total_dirac[0]);
            }
            if (needModel("bw")) {
                mt_tree->Branch("piminus_bw", &mt_data["piminus"].bw[0]);
                mt_tree->Branch("piminus_total_bw", &mt_data["piminus"].total_bw[0]);
            }
        }
        if (needParticle("pi0")) {
            mt_tree->Branch("pi0_prim", &mt_data["pi0"].prim[0]);
            if (needModel("dirac")) {
                mt_tree->Branch("pi0_dirac", &mt_data["pi0"].dirac[0]);
                mt_tree->Branch("pi0_total_dirac", &mt_data["pi0"].total_dirac[0]);
            }
            if (needModel("bw")) {
                mt_tree->Branch("pi0_bw", &mt_data["pi0"].bw[0]);
                mt_tree->Branch("pi0_total_bw", &mt_data["pi0"].total_bw[0]);
            }
        }

        // Wypełnienie drzewa – każdy wiersz odpowiada jednemu mt
        for (int i = 0; i < npoints_mt; ++i) {
            mt = mt_vals[i];
            mt_minus_mp = mt - mp;
            mt_minus_mn = mt - mn;
            mt_minus_mpi = mt - mpi;
            mt_minus_mpi0 = mt - mpi0;
            mt_tree->Fill();
        }
    }

    // -------------------------------------------------------------------------
    // Drzewo dla rozkładów rapidity (y>=0)
    // -------------------------------------------------------------------------
    TTree* y_tree = nullptr;
    if (need("rapidity")) {
        y_tree = new TTree("y_tree", "Rapidity distributions (y >= 0)");
        double y;
        y_tree->Branch("y", &y, "y/D");

        if (needParticle("proton")) {
            y_tree->Branch("proton_prim", &y_data["proton"].prim[0]);
            if (needModel("dirac")) {
                y_tree->Branch("proton_dirac", &y_data["proton"].dirac[0]);
                y_tree->Branch("proton_total_dirac", &y_data["proton"].total_dirac[0]);
            }
            if (needModel("bw")) {
                y_tree->Branch("proton_bw", &y_data["proton"].bw[0]);
                y_tree->Branch("proton_total_bw", &y_data["proton"].total_bw[0]);
            }
        }
        if (needParticle("neutron")) {
            y_tree->Branch("neutron_prim", &y_data["neutron"].prim[0]);
            if (needModel("dirac")) {
                y_tree->Branch("neutron_dirac", &y_data["neutron"].dirac[0]);
                y_tree->Branch("neutron_total_dirac", &y_data["neutron"].total_dirac[0]);
            }
            if (needModel("bw")) {
                y_tree->Branch("neutron_bw", &y_data["neutron"].bw[0]);
                y_tree->Branch("neutron_total_bw", &y_data["neutron"].total_bw[0]);
            }
        }
        if (needParticle("piplus")) {
            y_tree->Branch("piplus_prim", &y_data["piplus"].prim[0]);
            if (needModel("dirac")) {
                y_tree->Branch("piplus_dirac", &y_data["piplus"].dirac[0]);
                y_tree->Branch("piplus_total_dirac", &y_data["piplus"].total_dirac[0]);
            }
            if (needModel("bw")) {
                y_tree->Branch("piplus_bw", &y_data["piplus"].bw[0]);
                y_tree->Branch("piplus_total_bw", &y_data["piplus"].total_bw[0]);
            }
        }
        if (needParticle("piminus")) {
            y_tree->Branch("piminus_prim", &y_data["piminus"].prim[0]);
            if (needModel("dirac")) {
                y_tree->Branch("piminus_dirac", &y_data["piminus"].dirac[0]);
                y_tree->Branch("piminus_total_dirac", &y_data["piminus"].total_dirac[0]);
            }
            if (needModel("bw")) {
                y_tree->Branch("piminus_bw", &y_data["piminus"].bw[0]);
                y_tree->Branch("piminus_total_bw", &y_data["piminus"].total_bw[0]);
            }
        }
        if (needParticle("pi0")) {
            y_tree->Branch("pi0_prim", &y_data["pi0"].prim[0]);
            if (needModel("dirac")) {
                y_tree->Branch("pi0_dirac", &y_data["pi0"].dirac[0]);
                y_tree->Branch("pi0_total_dirac", &y_data["pi0"].total_dirac[0]);
            }
            if (needModel("bw")) {
                y_tree->Branch("pi0_bw", &y_data["pi0"].bw[0]);
                y_tree->Branch("pi0_total_bw", &y_data["pi0"].total_bw[0]);
            }
        }

        for (int i = 0; i < npoints_y_half; ++i) {
            y = y_half[i];
            y_tree->Fill();
        }
    }

    // -------------------------------------------------------------------------
    // Zapisz parametry
    // -------------------------------------------------------------------------
    TParameter<double>* p_scale = new TParameter<double>("scale_factor", scale_factor);
    TParameter<double>* p_Nexp = new TParameter<double>("N_proton_exp", N_proton_exp);
    p_scale->Write();
    p_Nexp->Write();

    outfile->Write();
    outfile->Close();

    timer.Stop();
    std::cout << "\n==========================================" << std::endl;
    std::cout << "Computation finished." << std::endl;
    std::cout << "Results saved to " << sel.output << std::endl;
    std::cout << "Time elapsed: " << timer.RealTime() << " seconds" << std::endl;
    std::cout << "==========================================" << std::endl;

    return 0;
}