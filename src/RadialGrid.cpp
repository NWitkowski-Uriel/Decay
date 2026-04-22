// src/RadialGrid.cpp
#include "RadialGrid.h"
#include "Constants.h"
#include <cmath>

RadialGrid& RadialGrid::instance() {
    static RadialGrid inst;
    return inst;
}

void RadialGrid::setNumPoints(int n) {
    if (n == n_) return;
    n_ = n;
    compute();
}

void RadialGrid::compute() {
    if (n_ <= 1) return;
    dr_ = R / (n_ - 1);
    r_.resize(n_);
    cosh_.resize(n_);
    sinh_.resize(n_);
    r2_.resize(n_);

    for (int i = 0; i < n_; ++i) {
        double r = i * dr_;
        r_[i] = r;
        cosh_[i] = std::cosh(H * r);
        sinh_[i] = std::sinh(H * r);
        r2_[i] = r * r;
    }
}