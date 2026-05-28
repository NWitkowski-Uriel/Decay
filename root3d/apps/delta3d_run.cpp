#include "Delta3D/Physics.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* exe) {
  std::cerr << "Usage: " << exe << " [--output DIR] [--serial] [--fast] [--production]\n\n"
            << "  --output DIR    output directory, default output/root3d\n"
            << "  --serial        disable OpenMP parallel loops\n"
            << "  --fast          reduce grids for quick smoke tests\n"
            << "  --production    increase grids for final tables and plots\n";
}

} // namespace

int main(int argc, char** argv) {
  auto run = delta3d::DefaultRunConfig();

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--output" && i + 1 < argc) {
      run.output_dir = argv[++i];
    } else if (arg == "--serial") {
      for (auto& cfg : run.parameter_sets) cfg.parallel = false;
    } else if (arg == "--fast") {
      for (auto& cfg : run.parameter_sets) {
        cfg.mass_grid = 80;
        cfg.p_grid = 80;
        cfg.r_grid = 60;
        cfg.mt_grid = 80;
        cfg.y_grid = 81;
      }
    } else if (arg == "--production") {
      for (auto& cfg : run.parameter_sets) {
        cfg.mass_grid = 600;
        cfg.p_grid = 320;
        cfg.r_grid = 220;
        cfg.mt_grid = 260;
        cfg.y_grid = 241;
      }
    } else if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      return EXIT_SUCCESS;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  try {
    std::cout << "Running Delta ROOT 3D analysis...\n";
    const auto result = delta3d::RunAnalysis(run);
    delta3d::WriteCsvTables(result, run.output_dir);
    delta3d::WriteRootFileAndPlots(result, run.output_dir);
    std::cout << "Done. Results written to " << run.output_dir << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& ex) {
    std::cerr << "ERROR: " << ex.what() << "\n";
    return EXIT_FAILURE;
  }
}
