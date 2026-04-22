// apps/compute_total.cpp
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TParameter.h"

#include "Constants.h"

bool hasNonZeroContent(TH1* h) {
    if (!h) return false;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        if (std::abs(h->GetBinContent(i)) > 1e-20) return true;
    }
    return false;
}

int main() {
    TFile f_prim("primordial.root", "READ");
    if (!f_prim.IsOpen()) {
        std::cerr << "Error: cannot open primordial.root\n";
        return 1;
    }

    TFile f_delta("delta.root", "READ");
    bool has_delta = f_delta.IsOpen();

    TFile f_out("total.root", "RECREATE");
    if (!f_out.IsOpen()) {
        std::cerr << "Error: cannot create total.root\n";
        return 1;
    }

    std::vector<std::string> particles = {"proton", "neutron", "piplus", "piminus", "pi0"};
    double total_proton_prim = 0.0;
    double total_proton_dirac = 0.0;

    for (const auto& p : particles) {
        // ----- mt histogram -----
        TH1D *h_prim_mt = (TH1D*)f_prim.Get(("h_prim_" + p + "_mt").c_str());
        if (!h_prim_mt) {
            std::cerr << "Warning: no primordial mt histogram for " << p << "\n";
            continue;
        }

        // Zawsze zapisz sumę pierwotną (total_prim)
        TH1D *h_total_prim_mt = (TH1D*)h_prim_mt->Clone(("h_total_prim_" + p + "_mt").c_str());
        h_total_prim_mt->Write();

        // Sprawdź, czy są wkłady z rozpadów
        TH1D *h_dirac_mt = nullptr, *h_bw_mt = nullptr;
        if (has_delta) {
            h_dirac_mt = (TH1D*)f_delta.Get(("h_delta_dirac_" + p + "_mt").c_str());
            h_bw_mt    = (TH1D*)f_delta.Get(("h_delta_bw_" + p + "_mt").c_str());
        }

        // Twórz total_dirac tylko jeśli istnieją dane Dirac i nie są puste
        if (h_dirac_mt && hasNonZeroContent(h_dirac_mt)) {
            TH1D *h_total_dirac_mt = (TH1D*)h_prim_mt->Clone(("h_total_dirac_" + p + "_mt").c_str());
            h_total_dirac_mt->Add(h_dirac_mt);
            h_total_dirac_mt->Write();
            delete h_total_dirac_mt;
        }

        // Podobnie dla BW
        if (h_bw_mt && hasNonZeroContent(h_bw_mt)) {
            TH1D *h_total_bw_mt = (TH1D*)h_prim_mt->Clone(("h_total_bw_" + p + "_mt").c_str());
            h_total_bw_mt->Add(h_bw_mt);
            h_total_bw_mt->Write();
            delete h_total_bw_mt;
        }

        delete h_total_prim_mt;

        // ----- y histogram -----
        TH1D *h_prim_y = (TH1D*)f_prim.Get(("h_prim_" + p + "_y").c_str());
        if (!h_prim_y) continue;

        TH1D *h_total_prim_y = (TH1D*)h_prim_y->Clone(("h_total_prim_" + p + "_y").c_str());
        h_total_prim_y->Write();

        if (has_delta) {
            TH1D *h_dirac_y = (TH1D*)f_delta.Get(("h_delta_dirac_" + p + "_y").c_str());
            if (h_dirac_y && hasNonZeroContent(h_dirac_y)) {
                TH1D *h_total_dirac_y = (TH1D*)h_prim_y->Clone(("h_total_dirac_" + p + "_y").c_str());
                h_total_dirac_y->Add(h_dirac_y);
                h_total_dirac_y->Write();
                if (p == "proton") {
                    for (int i = 1; i <= h_dirac_y->GetNbinsX(); ++i)
                        total_proton_dirac += h_dirac_y->GetBinContent(i) * h_dirac_y->GetBinWidth(i);
                }
                delete h_total_dirac_y;
            }

            TH1D *h_bw_y = (TH1D*)f_delta.Get(("h_delta_bw_" + p + "_y").c_str());
            if (h_bw_y && hasNonZeroContent(h_bw_y)) {
                TH1D *h_total_bw_y = (TH1D*)h_prim_y->Clone(("h_total_bw_" + p + "_y").c_str());
                h_total_bw_y->Add(h_bw_y);
                h_total_bw_y->Write();
                delete h_total_bw_y;
            }
        }

        // Całkowity yield protonu (pierwotny)
        if (p == "proton") {
            for (int i = 1; i <= h_prim_y->GetNbinsX(); ++i)
                total_proton_prim += h_prim_y->GetBinContent(i) * h_prim_y->GetBinWidth(i);
        }

        delete h_total_prim_y;
    }

    double total_proton_total_dirac = total_proton_prim + total_proton_dirac;
    double scale_factor = (total_proton_total_dirac > 0) ? N_proton_exp / total_proton_total_dirac : 1.0;

    std::cout << "Total proton (primordial): " << total_proton_prim << "\n";
    std::cout << "Total proton (Dirac decays): " << total_proton_dirac << "\n";
    std::cout << "Scale factor = " << scale_factor << "\n";

    TParameter<double> *p_scale = new TParameter<double>("scale_factor", scale_factor);
    p_scale->Write();

    f_out.Close();
    f_prim.Close();
    if (has_delta) f_delta.Close();

    return 0;
}