#include "Delta3D/PairFit.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace delta3d {
namespace {

constexpr double kMp = 938.2720813;
constexpr double kPiCharged = 139.57039;
constexpr double kPi = 3.141592653589793238462643383279502884;

std::vector<double> parse_numbers(const std::string& line) {
  std::string clean = line;
  for (char& c : clean) {
    if (c == ',' || c == ';' || c == '\t') c = ' ';
  }
  std::istringstream in(clean);
  std::vector<double> values;
  double x = 0.0;
  while (in >> x) values.push_back(x);
  return values;
}

double bw_shape_raw(double M, const ModelConfig& cfg, double mp, double mpi) {
  const double q = PStar(M, mp, mpi);
  const double q0 = PStar(cfg.m_delta, mp, mpi);
  if (q <= 0.0 || q0 <= 0.0 || M <= 0.0) return 0.0;
  const double gamma = cfg.gamma_delta * std::pow(q / q0, 3.0) * cfg.m_delta / M;
  const double den = (M*M - cfg.m_delta*cfg.m_delta) * (M*M - cfg.m_delta*cfg.m_delta)
                   + cfg.m_delta*cfg.m_delta*gamma*gamma;
  return den > 0.0 ? 2.0 * M * cfg.m_delta * gamma / den : 0.0;
}

double cugnon_shape_raw(double M, const ModelConfig& cfg, double mp, double mpi) {
  const double q = PStar(M, mp, mpi);
  const double q0 = PStar(cfg.m_delta, mp, mpi);
  if (q <= 0.0 || q0 <= 0.0 || M <= 0.0) return 0.0;
  const double R = std::max(0.0, cfg.cugnon_radius);
  const double barrier = (R > 0.0) ? ((1.0 + R*R*q0*q0) / (1.0 + R*R*q*q)) : 1.0;
  const double gamma = cfg.gamma_delta * std::pow(q / q0, 3.0) * cfg.m_delta / M * barrier;
  const double den = (M*M - cfg.m_delta*cfg.m_delta) * (M*M - cfg.m_delta*cfg.m_delta)
                   + cfg.m_delta*cfg.m_delta*gamma*gamma;
  return den > 0.0 ? 2.0 * M * cfg.m_delta * gamma / den : 0.0;
}

double phase_shift(double M, const ModelConfig& cfg, double mp, double mpi) {
  const double q = PStar(M, mp, mpi);
  if (q <= 0.0 || M <= 0.0) return 0.0;
  const double pi_ps = cfg.c1_ps*q*q + cfg.c2_ps*q*q*q*q;
  const double den = 3.0 * M * (M*M - cfg.m_delta*cfg.m_delta) * (1.0 + pi_ps);
  if (std::abs(den) < 1e-14) return kPi / 2.0;
  const double arg = -(2.0 * cfg.alpha0_ps * q*q*q) / den;
  return std::atan(arg);
}

double ps_shape_raw(double M, const ModelConfig& cfg, double mp, double mpi) {
  const double h = std::max(1.0e-3, 1.0e-5 * M);
  const double a = std::max(cfg.mass_min, M - h);
  const double b = std::min(cfg.mass_max, M + h);
  if (!(b > a)) return 0.0;
  const double rho = (phase_shift(b, cfg, mp, mpi) - phase_shift(a, cfg, mp, mpi)) / ((b - a) * kPi);
  return std::isfinite(rho) && rho > 0.0 ? rho : 0.0;
}

TF1 make_fit_function(MassModel model, const PairFitOptions& opt) {
  const double mp = kMp;
  const double mpi = kPiCharged;
  const std::string name = std::string("fit_") + ToString(model);
  TF1 f(name.c_str(), [=](double* x, double* p) {
    ModelConfig cfg = opt.seed;
    cfg.m_delta = p[1];
    cfg.gamma_delta = std::max(1.0, p[2]);
    cfg.cugnon_radius = std::max(0.0, p[5]);
    const double M = x[0];
    const double signal = PairLineShape(model, M, cfg, mp, mpi);
    const double background = p[3] + p[4] * (M - cfg.m_delta);
    return p[0] * signal + background;
  }, opt.fit_min, opt.fit_max, 6);

  f.SetParName(0, "A");
  f.SetParName(1, "M0");
  f.SetParName(2, "Gamma");
  f.SetParName(3, "B0");
  f.SetParName(4, "B1");
  f.SetParName(5, "Rcugnon");
  f.SetParameters(1.0, opt.seed.m_delta, opt.seed.gamma_delta, 0.0, 0.0, opt.seed.cugnon_radius);
  f.SetParLimits(1, opt.fit_min, opt.fit_max);
  f.SetParLimits(2, 1.0, 400.0);
  if (model == MassModel::Cugnon) {
    f.SetParLimits(5, 0.0, 0.05);
  } else {
    f.FixParameter(5, opt.seed.cugnon_radius);
  }
  return f;
}

void draw_pair_fit(const PairFitOptions& opt, const std::vector<PairFitResult>& results, const std::vector<PairDataPoint>& data) {
  fs::create_directories(opt.output_dir);
  TCanvas c("pairfit", "p pi+ mass spectrum fit", 1100, 800);
  TGraphErrors gr(static_cast<int>(data.size()));
  for (int i = 0; i < static_cast<int>(data.size()); ++i) {
    gr.SetPoint(i, data[i].mass, data[i].yield);
    gr.SetPointError(i, 0.0, std::max(data[i].error, 1.0e-12));
  }
  gr.SetTitle("p#pi^{+} pair mass spectrum;M_{p#pi^{+}} [MeV];dN/dM");
  gr.SetMarkerStyle(20);
  gr.Draw("AP");

  TLegend leg(0.60, 0.65, 0.88, 0.88);
  leg.AddEntry(&gr, "data", "p");
  int color = 2;
  for (const auto& r : results) {
    ModelConfig cfg = opt.seed;
    cfg.m_delta = r.mass;
    cfg.gamma_delta = r.gamma;
    auto* f = new TF1(("curve_" + r.model_name).c_str(), [=](double* x, double*) {
      return r.amplitude * PairLineShape(r.model, x[0], cfg, kMp, kPiCharged)
           + r.background0 + r.background1 * (x[0] - r.mass);
    }, opt.fit_min, opt.fit_max, 0);
    f->SetLineColor(color++);
    f->SetLineWidth(2);
    f->Draw("SAME");
    leg.AddEntry(f, r.model_name.c_str(), "l");
  }
  leg.Draw();
  c.SaveAs((opt.output_dir + "/pairfit_models.png").c_str());
}

} // namespace

std::vector<PairDataPoint> ReadPairSpectrum(const PairFitOptions& opt) {
  std::ifstream in(opt.input_path);
  if (!in) throw std::runtime_error("Cannot open pair spectrum file: " + opt.input_path);

  std::vector<PairDataPoint> data;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    const auto values = parse_numbers(line);
    const int required = opt.has_error_column ? std::max({opt.mass_column, opt.yield_column, opt.error_column}) : std::max(opt.mass_column, opt.yield_column);
    if (static_cast<int>(values.size()) <= required) continue;
    PairDataPoint p;
    p.mass = values[opt.mass_column] * (opt.input_mass_in_gev ? 1000.0 : 1.0);
    p.yield = values[opt.yield_column];
    p.error = opt.has_error_column ? std::abs(values[opt.error_column]) : std::sqrt(std::max(1.0, std::abs(p.yield)));
    if (std::isfinite(p.mass) && std::isfinite(p.yield) && p.mass > 0.0) data.push_back(p);
  }
  if (data.empty()) throw std::runtime_error("No valid pair spectrum rows were read from: " + opt.input_path);
  return data;
}

double PairLineShape(MassModel model, double M, const ModelConfig& cfg, double mp, double mpi) {
  switch (model) {
    case MassModel::Dirac:
      return 0.0;
    case MassModel::BreitWigner:
      return bw_shape_raw(M, cfg, mp, mpi);
    case MassModel::Cugnon:
      return cugnon_shape_raw(M, cfg, mp, mpi);
    case MassModel::PhaseShift:
      return ps_shape_raw(M, cfg, mp, mpi);
  }
  return 0.0;
}

std::vector<PairFitResult> FitPairSpectrum(const PairFitOptions& opt) {
  const auto data = ReadPairSpectrum(opt);
  TGraphErrors gr(static_cast<int>(data.size()));
  for (int i = 0; i < static_cast<int>(data.size()); ++i) {
    gr.SetPoint(i, data[i].mass, data[i].yield);
    gr.SetPointError(i, 0.0, std::max(data[i].error, 1.0e-12));
  }

  std::vector<PairFitResult> results;
  for (const auto model : {MassModel::BreitWigner, MassModel::Cugnon, MassModel::PhaseShift}) {
    TF1 f = make_fit_function(model, opt);
    gr.Fit(&f, "RQ0");
    PairFitResult r;
    r.model = model;
    r.model_name = ToString(model);
    r.amplitude = f.GetParameter(0);
    r.mass = f.GetParameter(1);
    r.gamma = f.GetParameter(2);
    r.background0 = f.GetParameter(3);
    r.background1 = f.GetParameter(4);
    r.chi2 = f.GetChisquare();
    r.ndf = f.GetNDF();
    r.chi2_ndf = r.ndf > 0 ? r.chi2 / r.ndf : std::numeric_limits<double>::infinity();
    results.push_back(r);
  }
  return results;
}

ModelConfig FittedConfigFromBestResult(const ModelConfig& seed, const std::vector<PairFitResult>& results) {
  if (results.empty()) return seed;
  const auto best = std::min_element(results.begin(), results.end(), [](const auto& a, const auto& b) {
    return a.chi2_ndf < b.chi2_ndf;
  });
  ModelConfig cfg = seed;
  cfg.tag = "fitted_" + best->model_name;
  cfg.m_delta = best->mass;
  cfg.gamma_delta = best->gamma;
  return cfg;
}

void WritePairFitOutputs(const PairFitOptions& opt, const std::vector<PairFitResult>& results, const std::vector<PairDataPoint>& data) {
  fs::create_directories(opt.output_dir);
  {
    std::ofstream out(opt.output_dir + "/pairfit_results.csv");
    out << "model,amplitude,mass_MeV,gamma_MeV,background0,background1,chi2,ndf,chi2_ndf\n";
    out << std::setprecision(12);
    for (const auto& r : results) {
      out << r.model_name << ',' << r.amplitude << ',' << r.mass << ',' << r.gamma << ','
          << r.background0 << ',' << r.background1 << ',' << r.chi2 << ',' << r.ndf << ',' << r.chi2_ndf << '\n';
    }
  }
  {
    const auto fitted = FittedConfigFromBestResult(opt.seed, results);
    std::ofstream out(opt.output_dir + "/fitted_config.txt");
    out << "# Paste these values into DefaultRunConfig() or load them in a later config reader.\n";
    out << std::setprecision(12);
    out << "tag " << fitted.tag << '\n';
    out << "m_delta " << fitted.m_delta << '\n';
    out << "gamma_delta " << fitted.gamma_delta << '\n';
    out << "mass_min " << fitted.mass_min << '\n';
    out << "mass_max " << fitted.mass_max << '\n';
  }
  draw_pair_fit(opt, results, data);
}

} // namespace delta3d
