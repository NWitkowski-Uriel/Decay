#pragma once

#include <string>
#include <vector>

namespace delta3d {

enum class MassModel { Dirac, BreitWigner, Cugnon, PhaseShift };

enum class ParticleKind { Proton, Neutron, PiPlus, PiZero, PiMinus, DeltaPP, DeltaP, Delta0, DeltaM };

struct Particle {
  ParticleKind kind{};
  std::string name;
  double mass = 0.0;      // MeV
  double degeneracy = 1.0;
  double chemical_mu = 0.0; // MeV
};

struct ModelConfig {
  std::string tag = "vacuum";
  double temperature = 49.0;      // MeV; match your Wolfram notebook before production runs
  double radial_flow_H = 0.0;     // 1/fm-like blast-wave parameter used in notebook kernels
  double radius = 8.0;            // fm; only relative normalization if parameters are not fitted
  double m_delta = 1232.0;        // MeV
  double gamma_delta = 117.0;     // MeV
  double mass_min = 938.2720813 + 139.57039;
  double mass_max = 1500.0;

  // Cugnon-inspired p-wave line-shape parameter used by PairFit.
  // The full 3D daughter-kernel propagation of the Cugnon model should be wired
  // together with the exact q-k-r correction kernel.
  double cugnon_radius = 0.0;

  // Pok-Man-Lo phase-shift parameters. Keep synchronized with the Wolfram PairFit notebook.
  double alpha0_ps = 0.2010;
  double c1_ps = 0.5320;
  double c2_ps = 0.0000;

  int mass_grid = 240;
  int p_grid = 180;
  int r_grid = 120;
  int mt_grid = 180;
  int y_grid = 161;
  double p_max_factor = 5.0;
  double mt_max_offset = 2000.0;  // MeV
  double y_max = 2.0;
  bool parallel = true;
};

struct RunConfig {
  std::string output_dir = "output/root3d";
  // The main 3D pipeline currently runs models with implemented propagation in Physics.cpp.
  // Cugnon is implemented in PairFit and will be enabled here after the exact q-k-r
  // daughter kernel is ported into the 3D correction spectrum.
  std::vector<MassModel> mass_models{MassModel::Dirac, MassModel::BreitWigner, MassModel::PhaseShift};
  std::vector<ModelConfig> parameter_sets{};
};

inline const char* ToString(MassModel m) {
  switch (m) {
    case MassModel::Dirac: return "Dirac";
    case MassModel::BreitWigner: return "BW";
    case MassModel::Cugnon: return "Cugnon";
    case MassModel::PhaseShift: return "PS";
  }
  return "Unknown";
}

} // namespace delta3d
