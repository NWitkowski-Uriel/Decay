#pragma once

#include "Delta3D/Config.h"
#include "Delta3D/Physics.h"

#include <string>
#include <vector>

namespace delta3d {

struct PairDataPoint {
  double mass = 0.0;   // MeV
  double yield = 0.0;
  double error = 1.0;
};

struct PairFitOptions {
  std::string input_path;
  std::string output_dir = "output/root3d/pairfit";
  double fit_min = 1100.0;
  double fit_max = 1400.0;
  double norm_min = 1077.8424713;
  double norm_max = 1500.0;
  bool input_mass_in_gev = true;
  bool has_error_column = false;
  int mass_column = 0;
  int yield_column = 1;
  int error_column = 2;
  ModelConfig seed{};
};

struct PairFitResult {
  MassModel model = MassModel::BreitWigner;
  std::string model_name;
  double amplitude = 1.0;
  double mass = 1232.0;
  double gamma = 117.0;
  double background0 = 0.0;
  double background1 = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
  double chi2_ndf = 0.0;
};

std::vector<PairDataPoint> ReadPairSpectrum(const PairFitOptions& opt);
double PairLineShape(MassModel model, double M, const ModelConfig& cfg, double mp, double mpi);
std::vector<PairFitResult> FitPairSpectrum(const PairFitOptions& opt);
ModelConfig FittedConfigFromBestResult(const ModelConfig& seed, const std::vector<PairFitResult>& results);
void WritePairFitOutputs(const PairFitOptions& opt, const std::vector<PairFitResult>& results, const std::vector<PairDataPoint>& data);

} // namespace delta3d
