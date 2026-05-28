#pragma once

#include "Delta3D/Config.h"

#include <functional>
#include <map>
#include <string>
#include <vector>

class TDirectory;

namespace delta3d {

struct SpectrumPoint {
  double x = 0.0;
  double value = 0.0;
};

struct YieldRow {
  std::string parameter_tag;
  std::string mass_model;
  std::string particle;
  double primordial = 0.0;
  double correction = 0.0;
  double total = 0.0;
};

struct BranchTestRow {
  std::string parameter_tag;
  std::string mass_model;
  std::string channel;
  std::string observable;
  double expected = 1.0;
  double value = 0.0;
};

struct ParticleSpectra {
  std::string particle;
  std::vector<SpectrumPoint> mt_primordial;
  std::vector<SpectrumPoint> mt_correction;
  std::vector<SpectrumPoint> mt_total;
  std::vector<SpectrumPoint> y_primordial;
  std::vector<SpectrumPoint> y_correction;
  std::vector<SpectrumPoint> y_total;
};

struct RunResult {
  std::vector<YieldRow> yields;
  std::vector<BranchTestRow> tests;
  std::map<std::string, ParticleSpectra> spectra;
};

double PStar(double M, double m1, double m2);
double KMinus(double q, double M, double m1, double m2);
double KPlus(double q, double M, double m1, double m2);

double Integrate1D(const std::function<double(double)>& f, double a, double b, int n);
double Integrate2D(const std::function<double(double,double)>& f, double ax, double bx, double ay, double by, int nx, int ny);
double Integrate3D(const std::function<double(double,double,double)>& f, double ax, double bx, double ay, double by, double az, double bz, int nx, int ny, int nz);

double SpectralWeight(MassModel model, double M, const ModelConfig& cfg, double m1, double m2);
double ThermalKernel3D(double p, double r, double mass, double mu, double g, const ModelConfig& cfg);
double ThermalYieldDirac(const Particle& particle, const ModelConfig& cfg);
double ThermalYieldMassDistributed(const Particle& delta, MassModel model, const ModelConfig& cfg, double daughter1, double daughter2);
double PrimordialMt(double mt, double y, const Particle& particle, const ModelConfig& cfg);
double PrimordialRapidity(double y, const Particle& particle, const ModelConfig& cfg);

RunConfig DefaultRunConfig();
RunResult RunAnalysis(const RunConfig& run);
void WriteCsvTables(const RunResult& result, const std::string& output_dir);
void WriteRootFileAndPlots(const RunResult& result, const std::string& output_dir);

} // namespace delta3d
