// include/RadialGrid.h
#ifndef RADIALGRID_H
#define RADIALGRID_H

#include <vector>

class RadialGrid {
public:
    // Singleton pattern – one instance per program
    static RadialGrid& instance();

    // Number of radial points (can be set once)
    void setNumPoints(int n);

    // Accessors
    const std::vector<double>& r()   const { return r_; }
    const std::vector<double>& cosh() const { return cosh_; }
    const std::vector<double>& sinh() const { return sinh_; }
    const std::vector<double>& r2()  const { return r2_; }
    double dr() const { return dr_; }

private:
    RadialGrid() : n_(0), dr_(0.0) {}
    void compute();

    int n_;
    double dr_;
    std::vector<double> r_;
    std::vector<double> cosh_;
    std::vector<double> sinh_;
    std::vector<double> r2_;
};

#endif // RADIALGRID_H