// src/Plotting.cpp
#include "Plotting.h"
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <iostream>
#include <filesystem>

void DrawComparisonPlot(
    const std::string& filename,
    const std::string& title,
    const std::string& x_title,
    const std::string& y_title,
    const std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>>& series,
    const std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>>& exp_data,
    bool logy,
    double x_min, double x_max,
    double y_min, double y_max)
{
    std::filesystem::path path(filename);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    bool has_series = false;
    for (const auto& s : series) {
        auto x = std::get<0>(s);
        auto y = std::get<1>(s);
        if (x && y && !x->empty() && !y->empty()) {
            has_series = true;
            break;
        }
    }
    if (!has_series) {
        std::cerr << "Warning: No valid series data for plot " << filename << " – skipping." << std::endl;
        return;
    }

    TCanvas* c = new TCanvas(filename.c_str(), title.c_str(), 1200, 800);
    if (logy) c->SetLogy();
    c->SetGrid();
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(title.c_str());

    for (const auto& s : series) {
        auto x = std::get<0>(s);
        auto y = std::get<1>(s);
        int color = std::get<2>(s);
        int style = std::get<3>(s);
        std::string leg_title = std::get<4>(s);

        if (!x || !y || x->empty() || y->empty()) continue;
        TGraph* g = new TGraph(x->size(), x->data(), y->data());
        g->SetLineColor(color);
        g->SetLineWidth(2);
        g->SetLineStyle(style);
        g->SetTitle(leg_title.c_str());
        mg->Add(g);
    }

    mg->Draw("AL");
    mg->GetXaxis()->SetTitle(x_title.c_str());
    if (x_max > x_min) mg->GetXaxis()->SetLimits(x_min, x_max);
    mg->GetYaxis()->SetTitle(y_title.c_str());
    if (y_max > y_min) mg->GetYaxis()->SetRangeUser(y_min, y_max);

    c->SaveAs((filename + ".pdf").c_str());
    c->SaveAs((filename + ".png").c_str());

    delete c;
}
