#include <functional>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "TFile.h"
#include "TH1D.h"

#include "Constants.h"
#include "DecayFunctions.h"
#include "RadialGrid.h"
#include "Thermodynamics.h"

namespace {

// unchanged args and parsing...

void fill_mt_hist(TH1D& hist, double mass, const std::function<double(double)>& fn) {
    const int n = hist.GetNbinsX();
    std::vector<double> tmp(n, 0.0);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double mt = hist.GetBinCenter(i + 1);
        tmp[i] = (mt > mass) ? fn(mt) : 0.0;
    }

    for (int i = 0; i < n; ++i) hist.SetBinContent(i + 1, tmp[i]);
}

void fill_y_hist(TH1D& hist, const std::function<double(double)>& fn) {
    const int n = hist.GetNbinsX();
    const int half = n / 2;
    std::vector<double> tmp(n, 0.0);

#pragma omp parallel for
    for (int i = 0; i <= half; ++i) {
        double y = hist.GetBinCenter(i + 1);
        double v = fn(std::abs(y));
        tmp[i] = v;
        tmp[n - i - 1] = v;
    }

    for (int i = 0; i < n; ++i) hist.SetBinContent(i + 1, tmp[i]);
}

} // namespace

// rest unchanged
