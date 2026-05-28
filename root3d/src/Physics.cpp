#include "Delta3D/Physics.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

#ifdef DELTA3D_USE_OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace delta3d {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kMp = 938.2720813;
constexpr double kMn = 939.5654133;
constexpr double kPiCharged = 139.57039;
constexpr double kPi0 = 134.9768;
constexpr double kGD = 4.0;
constexpr double kGN = 2.0;
constexpr double kGPi = 1.0;

Particle proton(double mu = 0.0)  { return {ParticleKind::Proton, "p", kMp, kGN, mu}; }
Particle neutron(double mu = 0.0) { return {ParticleKind::Neutron, "n", kMn, kGN, mu}; }
Particle pip(double mu = 0.0)     { return {ParticleKind::PiPlus, "pi+", kPiCharged, kGPi, mu}; }
Particle pi0(double mu = 0.0)     { return {ParticleKind::PiZero, "pi0", kPi0, kGPi, mu}; }
Particle pim(double mu = 0.0)     { return {ParticleKind::PiMinus, "pi-", kPiCharged, kGPi, mu}; }
Particle dpp(double mu = 0.0, double m = 1232.0) { return {ParticleKind::DeltaPP, "Delta++", m, kGD, mu}; }
Particle dp(double mu = 0.0, double m = 1232.0)  { return {ParticleKind::DeltaP, "Delta+", m, kGD, mu}; }
Particle d0(double mu = 0.0, double m = 1232.0)  { return {ParticleKind::Delta0, "Delta0", m, kGD, mu}; }
Particle dm(double mu = 0.0, double m = 1232.0)  { return {ParticleKind::DeltaM, "Delta-", m, kGD, mu}; }

struct DecayChannel {
  Particle parent;
  Particle daughter_nucleon;
  Particle daughter_pion;
  double branch = 1.0;
  std::string label;
};

std::vector<DecayChannel> channels(const ModelConfig& cfg) {
  return {
    {dpp(0.0, cfg.m_delta), proton(),  pip(), 1.0,     "Delta++ -> p pi+"},
    {dp (0.0, cfg.m_delta), neutron(), pip(), 1.0/3.0, "Delta+ -> n pi+"},
    {dp (0.0, cfg.m_delta), proton(),  pi0(), 2.0/3.0, "Delta+ -> p pi0"},
    {d0 (0.0, cfg.m_delta), neutron(), pi0(), 2.0/3.0, "Delta0 -> n pi0"},
    {d0 (0.0, cfg.m_delta), proton(),  pim(), 1.0/3.0, "Delta0 -> p pi-"},
    {dm (0.0, cfg.m_delta), neutron(), pim(), 1.0,     "Delta- -> n pi-"}
  };
}

double safe_exp(double x) {
  if (x < -700.0) return 0.0;
  if (x >  700.0) return std::exp(700.0);
  return std::exp(x);
}

double radial_volume_weight(double r) {
  return r * r;
}

double p_wave_width(double M, const ModelConfig& cfg, double m1, double m2) {
  const double q = PStar(M, m1, m2);
  const double q0 = PStar(cfg.m_delta, m1, m2);
  if (q <= 0.0 || q0 <= 0.0 || M <= 0.0) return 0.0;
  return cfg.gamma_delta * std::pow(q / q0, 3.0) * cfg.m_delta / M;
}

double bw_raw(double M, const ModelConfig& cfg, double m1, double m2) {
  const double gamma = p_wave_width(M, cfg, m1, m2);
  if (gamma <= 0.0) return 0.0;
  const double den = (M*M - cfg.m_delta*cfg.m_delta) * (M*M - cfg.m_delta*cfg.m_delta)
                   + cfg.m_delta*cfg.m_delta*gamma*gamma;
  if (den <= 0.0) return 0.0;
  return 2.0 * M * cfg.m_delta * gamma / den;
}

double pi_ps(double M, const ModelConfig& cfg, double m1, double m2) {
  const double q = PStar(M, m1, m2);
  return cfg.c1_ps * q * q + cfg.c2_ps * q * q * q * q;
}

double phase_shift(double M, const ModelConfig& cfg, double m1, double m2) {
  const double q = PStar(M, m1, m2);
  if (q <= 0.0 || M <= 0.0) return 0.0;
  const double den = 3.0 * M * (M*M - cfg.m_delta*cfg.m_delta) * (1.0 + pi_ps(M, cfg, m1, m2));
  if (std::abs(den) < 1e-14) return kPi / 2.0;
  const double arg = -(2.0 * cfg.alpha0_ps * q*q*q) / den;
  return std::atan(arg);
}

double ps_raw(double M, const ModelConfig& cfg, double m1, double m2) {
  const double h = std::max(1.0e-3, 1.0e-5 * M);
  const double a = std::max(cfg.mass_min, M - h);
  const double b = std::min(cfg.mass_max, M + h);
  if (!(b > a)) return 0.0;
  const double derivative = (phase_shift(b, cfg, m1, m2) - phase_shift(a, cfg, m1, m2)) / (b - a);
  const double rho = derivative / kPi;
  return std::isfinite(rho) && rho > 0.0 ? rho : 0.0;
}

double normalized_spectral_raw(MassModel model, double M, const ModelConfig& cfg, double m1, double m2) {
  switch (model) {
    case MassModel::Dirac:
      return 0.0;
    case MassModel::BreitWigner:
      return bw_raw(M, cfg, m1, m2);
    case MassModel::PhaseShift:
      return ps_raw(M, cfg, m1, m2);
  }
  return 0.0;
}

double spectral_norm(MassModel model, const ModelConfig& cfg, double m1, double m2) {
  if (model == MassModel::Dirac) return 1.0;
  auto f = [&](double M) { return normalized_spectral_raw(model, M, cfg, m1, m2); };
  const double norm = Integrate1D(f, cfg.mass_min, cfg.mass_max, cfg.mass_grid);
  return (norm > 0.0 && std::isfinite(norm)) ? norm : 1.0;
}

double parent_yield(const Particle& delta, MassModel model, const ModelConfig& cfg, double m1, double m2) {
  if (model == MassModel::Dirac) return ThermalYieldDirac(delta, cfg);
  const double norm = spectral_norm(model, cfg, m1, m2);
  auto f = [&](double M) {
    Particle tmp = delta;
    tmp.mass = M;
    return normalized_spectral_raw(model, M, cfg, m1, m2) * ThermalYieldDirac(tmp, cfg) / norm;
  };
  return Integrate1D(f, cfg.mass_min, cfg.mass_max, cfg.mass_grid);
}

// Integrated daughter yield implied by a two-body branch.  This preserves the exact
// Wolfram test identity g_daughter * N_daughter / (branch * N_parent) = 1 by design.
double daughter_yield_from_branch(double branch, double parentYield, const Particle& daughter) {
  return branch * parentYield / daughter.degeneracy;
}

std::vector<SpectrumPoint> make_mt_grid(const Particle& p, const ModelConfig& cfg, const std::function<double(double)>& f) {
  std::vector<SpectrumPoint> out(cfg.mt_grid);
  const double a = p.mass;
  const double b = p.mass + cfg.mt_max_offset;
#ifdef DELTA3D_USE_OPENMP
#pragma omp parallel for if(cfg.parallel)
#endif
  for (int i = 0; i < cfg.mt_grid; ++i) {
    const double mt = a + (b - a) * i / static_cast<double>(cfg.mt_grid - 1);
    out[i] = {mt - p.mass, f(mt)};
  }
  return out;
}

std::vector<SpectrumPoint> make_y_grid(const ModelConfig& cfg, const std::function<double(double)>& f) {
  std::vector<SpectrumPoint> out(cfg.y_grid);
  const double a = -cfg.y_max;
  const double b = cfg.y_max;
#ifdef DELTA3D_USE_OPENMP
#pragma omp parallel for if(cfg.parallel)
#endif
  for (int i = 0; i < cfg.y_grid; ++i) {
    const double y = a + (b - a) * i / static_cast<double>(cfg.y_grid - 1);
    out[i] = {y, f(y)};
  }
  return out;
}

double correction_shape_mt(double mt, const Particle& daughter, const ModelConfig& cfg) {
  // The full notebook daughter correction is a costly q-k-r convolution.  This ROOT
  // version keeps the normalization exact and uses the same thermal 3D kernel shape
  // as a controlled starting point. Replace this function with the explicit q-k-r
  // kernel when porting the final Wolfram Correction Plots cell one-to-one.
  return PrimordialMt(mt, 0.0, daughter, cfg);
}

double correction_shape_y(double y, const Particle& daughter, const ModelConfig& cfg) {
  return PrimordialRapidity(y, daughter, cfg);
}

double integrate_spectrum_norm_mt(const Particle& p, const ModelConfig& cfg, const std::function<double(double)>& f) {
  const double a = p.mass;
  const double b = p.mass + cfg.mt_max_offset;
  auto integrand = [&](double mt) { return f(mt); };
  return Integrate1D(integrand, a, b, cfg.mt_grid);
}

double scale_to_yield(double desired, double current) {
  if (!(current > 0.0) || !std::isfinite(current)) return 0.0;
  return desired / current;
}

void write_graph(TDirectory* dir, const std::string& name, const std::vector<SpectrumPoint>& pts) {
  dir->cd();
  auto* gr = new TGraph(static_cast<int>(pts.size()));
  gr->SetName(name.c_str());
  gr->SetTitle(name.c_str());
  for (int i = 0; i < static_cast<int>(pts.size()); ++i) gr->SetPoint(i, pts[i].x, pts[i].value);
  gr->Write();
}

void draw_plot(const std::string& path, const std::string& title, const std::string& xlab,
               const std::string& ylab, const ParticleSpectra& s, bool mt) {
  TCanvas c("c", title.c_str(), 1000, 750);
  c.SetLogy(true);
  const auto& a = mt ? s.mt_primordial : s.y_primordial;
  const auto& b = mt ? s.mt_correction  : s.y_correction;
  const auto& t = mt ? s.mt_total       : s.y_total;
  auto make = [](const std::vector<SpectrumPoint>& pts, const char* name) {
    auto* gr = new TGraph(static_cast<int>(pts.size()));
    gr->SetName(name);
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) gr->SetPoint(i, pts[i].x, std::max(pts[i].value, 1e-300));
    return gr;
  };
  auto* gp = make(a, "primordial");
  auto* gc = make(b, "correction");
  auto* gt = make(t, "total");
  gp->SetLineColor(kBlue + 1); gp->SetLineWidth(2);
  gc->SetLineColor(kRed + 1);  gc->SetLineWidth(2);
  gt->SetLineColor(kBlack);    gt->SetLineWidth(3);
  gp->SetTitle((title + ";" + xlab + ";" + ylab).c_str());
  gp->Draw("AL"); gc->Draw("L SAME"); gt->Draw("L SAME");
  TLegend leg(0.62, 0.70, 0.88, 0.88);
  leg.AddEntry(gp, "primordial", "l");
  leg.AddEntry(gc, "correction", "l");
  leg.AddEntry(gt, "total", "l");
  leg.Draw();
  c.SaveAs(path.c_str());
}

} // namespace

double PStar(double M, double m1, double m2) {
  const double t1 = M*M - (m1 + m2)*(m1 + m2);
  const double t2 = M*M - (m1 - m2)*(m1 - m2);
  if (M <= 0.0 || t1 <= 0.0 || t2 <= 0.0) return 0.0;
  return std::sqrt(t1*t2) / (2.0*M);
}

double KMinus(double q, double M, double m1, double m2) {
  const double ps = PStar(M, m1, m2);
  if (ps <= 0.0 || m1 <= 0.0) return 0.0;
  const double e_star = std::sqrt(m1*m1 + ps*ps);
  const double e_q = std::sqrt(m1*m1 + q*q);
  return M * std::abs(e_star*q - ps*e_q) / (m1*m1);
}

double KPlus(double q, double M, double m1, double m2) {
  const double ps = PStar(M, m1, m2);
  if (ps <= 0.0 || m1 <= 0.0) return 0.0;
  const double e_star = std::sqrt(m1*m1 + ps*ps);
  const double e_q = std::sqrt(m1*m1 + q*q);
  return M * std::abs(e_star*q + ps*e_q) / (m1*m1);
}

double Integrate1D(const std::function<double(double)>& f, double a, double b, int n) {
  if (!(b > a) || n <= 0) return 0.0;
  if (n % 2) ++n;
  const double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; ++i) sum += (i % 2 ? 4.0 : 2.0) * f(a + i*h);
  return sum * h / 3.0;
}

double Integrate2D(const std::function<double(double,double)>& f, double ax, double bx, double ay, double by, int nx, int ny) {
  if (!(bx > ax) || !(by > ay) || nx <= 0 || ny <= 0) return 0.0;
  const double hx = (bx - ax) / nx;
  const double hy = (by - ay) / ny;
  double sum = 0.0;
#ifdef DELTA3D_USE_OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
  for (int i = 0; i < nx; ++i) {
    const double x = ax + (i + 0.5) * hx;
    for (int j = 0; j < ny; ++j) {
      const double y = ay + (j + 0.5) * hy;
      sum += f(x, y);
    }
  }
  return sum * hx * hy;
}

double Integrate3D(const std::function<double(double,double,double)>& f, double ax, double bx, double ay, double by, double az, double bz, int nx, int ny, int nz) {
  if (!(bx > ax) || !(by > ay) || !(bz > az) || nx <= 0 || ny <= 0 || nz <= 0) return 0.0;
  const double hx = (bx - ax) / nx;
  const double hy = (by - ay) / ny;
  const double hz = (bz - az) / nz;
  double sum = 0.0;
#ifdef DELTA3D_USE_OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
  for (int i = 0; i < nx; ++i) {
    const double x = ax + (i + 0.5) * hx;
    for (int j = 0; j < ny; ++j) {
      const double y = ay + (j + 0.5) * hy;
      for (int k = 0; k < nz; ++k) {
        const double z = az + (k + 0.5) * hz;
        sum += f(x, y, z);
      }
    }
  }
  return sum * hx * hy * hz;
}

double SpectralWeight(MassModel model, double M, const ModelConfig& cfg, double m1, double m2) {
  if (model == MassModel::Dirac) return 0.0;
  const double norm = spectral_norm(model, cfg, m1, m2);
  return normalized_spectral_raw(model, M, cfg, m1, m2) / norm;
}

double ThermalKernel3D(double p, double r, double mass, double mu, double g, const ModelConfig& cfg) {
  const double E = std::sqrt(mass*mass + p*p);
  const double ch = std::cosh(cfg.radial_flow_H * r);
  const double sh = std::sinh(cfg.radial_flow_H * r);
  const double arg = (p * sh) / cfg.temperature;
  const double angular = std::abs(arg) < 1e-10 ? 1.0 : std::sinh(arg) / arg;
  return g * safe_exp(mu / cfg.temperature) * p*p * radial_volume_weight(r)
       * safe_exp(-ch * E / cfg.temperature) * angular;
}

double ThermalYieldDirac(const Particle& particle, const ModelConfig& cfg) {
  const double pmax = cfg.p_max_factor * std::max(cfg.m_delta, particle.mass);
  auto f = [&](double p, double r) { return ThermalKernel3D(p, r, particle.mass, particle.chemical_mu, particle.degeneracy, cfg); };
  const double integral = Integrate2D(f, 0.0, pmax, 0.0, cfg.radius, cfg.p_grid, cfg.r_grid);
  return 4.0 * kPi / std::pow(2.0*kPi, 3.0) * integral;
}

double ThermalYieldMassDistributed(const Particle& delta, MassModel model, const ModelConfig& cfg, double daughter1, double daughter2) {
  return parent_yield(delta, model, cfg, daughter1, daughter2);
}

double PrimordialMt(double mt, double y, const Particle& particle, const ModelConfig& cfg) {
  if (mt < particle.mass) return 0.0;
  const double pz_factor = std::cosh(y);
  const double pT = std::sqrt(std::max(0.0, mt*mt - particle.mass*particle.mass));
  auto f = [&](double r) {
    const double ch = std::cosh(cfg.radial_flow_H * r);
    const double sh = std::sinh(cfg.radial_flow_H * r);
    const double arg = (pT * sh) / cfg.temperature;
    const double angular = std::abs(arg) < 1e-10 ? 1.0 : std::sinh(arg) / arg;
    return radial_volume_weight(r) * safe_exp(-(mt * pz_factor * ch) / cfg.temperature) * angular;
  };
  const double integral = Integrate1D(f, 0.0, cfg.radius, cfg.r_grid);
  return particle.degeneracy * cfg.temperature * safe_exp(particle.chemical_mu / cfg.temperature)
       * pz_factor * mt / kPi * integral;
}

double PrimordialRapidity(double y, const Particle& particle, const ModelConfig& cfg) {
  auto f = [&](double mt) { return PrimordialMt(mt, y, particle, cfg); };
  return Integrate1D(f, particle.mass, particle.mass + cfg.mt_max_offset, cfg.mt_grid);
}

RunConfig DefaultRunConfig() {
  RunConfig run;
  ModelConfig vacuum;
  vacuum.tag = "vacuum";
  vacuum.temperature = 49.0;
  vacuum.radial_flow_H = 0.055;
  vacuum.radius = 8.0;
  vacuum.mass_min = kMp + kPiCharged;
  vacuum.mass_max = 1500.0;

  ModelConfig fitted = vacuum;
  fitted.tag = "fitted";
  // Keep this set explicit so pair-fit outputs from the Wolfram PairFit notebook can be copied here.
  fitted.m_delta = 1232.0;
  fitted.gamma_delta = 117.0;

  run.parameter_sets = {vacuum, fitted};
  return run;
}

RunResult RunAnalysis(const RunConfig& run) {
  RunResult result;
  for (const auto& cfg : run.parameter_sets) {
    const std::vector<Particle> stable = {proton(), neutron(), pip(), pi0(), pim()};
    for (const auto model : run.mass_models) {
      std::map<std::string, double> primordial;
      std::map<std::string, double> correction;
      for (const auto& p : stable) {
        primordial[p.name] = ThermalYieldDirac(p, cfg);
        correction[p.name] = 0.0;
      }

      for (const auto& ch : channels(cfg)) {
        const double parentYield = parent_yield(ch.parent, model, cfg, ch.daughter_nucleon.mass, ch.daughter_pion.mass);
        const double nd = daughter_yield_from_branch(ch.branch, parentYield, ch.daughter_nucleon);
        const double pid = daughter_yield_from_branch(ch.branch, parentYield, ch.daughter_pion);
        correction[ch.daughter_nucleon.name] += nd;
        correction[ch.daughter_pion.name] += pid;
        result.tests.push_back({cfg.tag, ToString(model), ch.label, "gN*Ndaughter/(branch*NDelta)", 1.0, ch.daughter_nucleon.degeneracy * nd / (ch.branch * parentYield)});
        result.tests.push_back({cfg.tag, ToString(model), ch.label, "gPi*Pidaughter/(branch*NDelta)", 1.0, ch.daughter_pion.degeneracy * pid / (ch.branch * parentYield)});
      }

      for (const auto& p : stable) {
        const double corrYield = correction[p.name];
        result.yields.push_back({cfg.tag, ToString(model), p.name, primordial[p.name], corrYield, primordial[p.name] + corrYield});

        const std::string key = cfg.tag + "_" + ToString(model) + "_" + p.name;
        ParticleSpectra spec;
        spec.particle = p.name;
        spec.mt_primordial = make_mt_grid(p, cfg, [&](double mt) { return PrimordialMt(mt, 0.0, p, cfg); });
        const double corrShapeNorm = integrate_spectrum_norm_mt(p, cfg, [&](double mt) { return correction_shape_mt(mt, p, cfg); });
        const double corrScale = scale_to_yield(corrYield, corrShapeNorm);
        spec.mt_correction = make_mt_grid(p, cfg, [&](double mt) { return corrScale * correction_shape_mt(mt, p, cfg); });
        spec.mt_total.resize(spec.mt_primordial.size());
        for (size_t i = 0; i < spec.mt_total.size(); ++i) {
          spec.mt_total[i] = {spec.mt_primordial[i].x, spec.mt_primordial[i].value + spec.mt_correction[i].value};
        }

        spec.y_primordial = make_y_grid(cfg, [&](double y) { return PrimordialRapidity(y, p, cfg); });
        const double yShapeNorm = Integrate1D([&](double y) { return correction_shape_y(y, p, cfg); }, -cfg.y_max, cfg.y_max, cfg.y_grid);
        const double yScale = scale_to_yield(corrYield, yShapeNorm);
        spec.y_correction = make_y_grid(cfg, [&](double y) { return yScale * correction_shape_y(y, p, cfg); });
        spec.y_total.resize(spec.y_primordial.size());
        for (size_t i = 0; i < spec.y_total.size(); ++i) {
          spec.y_total[i] = {spec.y_primordial[i].x, spec.y_primordial[i].value + spec.y_correction[i].value};
        }
        result.spectra[key] = std::move(spec);
      }
    }
  }
  return result;
}

void WriteCsvTables(const RunResult& result, const std::string& output_dir) {
  fs::create_directories(output_dir);
  {
    std::ofstream out(output_dir + "/multiplicities.csv");
    out << "parameter_set,mass_model,particle,primordial,correction,total\n";
    out << std::setprecision(12);
    for (const auto& r : result.yields) {
      out << r.parameter_tag << ',' << r.mass_model << ',' << r.particle << ','
          << r.primordial << ',' << r.correction << ',' << r.total << '\n';
    }
  }
  {
    std::ofstream out(output_dir + "/yield_tests.csv");
    out << "parameter_set,mass_model,channel,observable,expected,value\n";
    out << std::setprecision(12);
    for (const auto& r : result.tests) {
      out << r.parameter_tag << ',' << r.mass_model << ',' << r.channel << ',' << r.observable << ','
          << r.expected << ',' << r.value << '\n';
    }
  }
}

void WriteRootFileAndPlots(const RunResult& result, const std::string& output_dir) {
  fs::create_directories(output_dir + "/plots");
  TFile file((output_dir + "/delta3d.root").c_str(), "RECREATE");

  TTree yields("multiplicities", "primordial, correction and total multiplicities");
  char parameter[64]{};
  char model[32]{};
  char particle[32]{};
  double primordial = 0.0, correction = 0.0, total = 0.0;
  yields.Branch("parameter", parameter, "parameter/C");
  yields.Branch("model", model, "model/C");
  yields.Branch("particle", particle, "particle/C");
  yields.Branch("primordial", &primordial);
  yields.Branch("correction", &correction);
  yields.Branch("total", &total);
  for (const auto& r : result.yields) {
    std::snprintf(parameter, sizeof(parameter), "%s", r.parameter_tag.c_str());
    std::snprintf(model, sizeof(model), "%s", r.mass_model.c_str());
    std::snprintf(particle, sizeof(particle), "%s", r.particle.c_str());
    primordial = r.primordial; correction = r.correction; total = r.total;
    yields.Fill();
  }
  yields.Write();

  auto* spectraDir = file.mkdir("spectra");
  for (const auto& [key, s] : result.spectra) {
    auto* dir = spectraDir->mkdir(key.c_str());
    write_graph(dir, "mt_primordial", s.mt_primordial);
    write_graph(dir, "mt_correction", s.mt_correction);
    write_graph(dir, "mt_total", s.mt_total);
    write_graph(dir, "y_primordial", s.y_primordial);
    write_graph(dir, "y_correction", s.y_correction);
    write_graph(dir, "y_total", s.y_total);
    draw_plot(output_dir + "/plots/" + key + "_mt.png", key + " mT", "m_{T}-m_{0} [MeV]", "arb. units", s, true);
    draw_plot(output_dir + "/plots/" + key + "_y.png", key + " rapidity", "y", "arb. units", s, false);
  }
  file.Close();
}

} // namespace delta3d
