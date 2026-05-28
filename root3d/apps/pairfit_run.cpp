#include "Delta3D/PairFit.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void usage(const char* exe) {
  std::cerr << "Usage: " << exe << " --input FILE [options]\n\n"
            << "Options:\n"
            << "  --output DIR        output directory, default output/root3d/pairfit\n"
            << "  --fit-min MeV       lower fit limit, default 1100\n"
            << "  --fit-max MeV       upper fit limit, default 1400\n"
            << "  --mass-mev          input mass column is already in MeV; default assumes GeV\n"
            << "  --has-errors        read error column\n"
            << "  --mass-col N        zero-based mass column, default 0\n"
            << "  --yield-col N       zero-based yield column, default 1\n"
            << "  --error-col N       zero-based error column, default 2\n";
}

} // namespace

int main(int argc, char** argv) {
  delta3d::PairFitOptions opt;
  opt.seed = delta3d::DefaultRunConfig().parameter_sets.front();

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--input" && i + 1 < argc) opt.input_path = argv[++i];
    else if (arg == "--output" && i + 1 < argc) opt.output_dir = argv[++i];
    else if (arg == "--fit-min" && i + 1 < argc) opt.fit_min = std::stod(argv[++i]);
    else if (arg == "--fit-max" && i + 1 < argc) opt.fit_max = std::stod(argv[++i]);
    else if (arg == "--mass-mev") opt.input_mass_in_gev = false;
    else if (arg == "--has-errors") opt.has_error_column = true;
    else if (arg == "--mass-col" && i + 1 < argc) opt.mass_column = std::stoi(argv[++i]);
    else if (arg == "--yield-col" && i + 1 < argc) opt.yield_column = std::stoi(argv[++i]);
    else if (arg == "--error-col" && i + 1 < argc) opt.error_column = std::stoi(argv[++i]);
    else if (arg == "--help" || arg == "-h") { usage(argv[0]); return EXIT_SUCCESS; }
    else { std::cerr << "Unknown or incomplete option: " << arg << "\n"; usage(argv[0]); return EXIT_FAILURE; }
  }

  if (opt.input_path.empty()) {
    std::cerr << "Missing required --input FILE\n";
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  try {
    const auto data = delta3d::ReadPairSpectrum(opt);
    const auto results = delta3d::FitPairSpectrum(opt);
    delta3d::WritePairFitOutputs(opt, results, data);
    std::cout << "Pair fit completed. Results written to " << opt.output_dir << "\n";
    for (const auto& r : results) {
      std::cout << r.model_name << ": M0=" << r.mass << " MeV, Gamma=" << r.gamma
                << " MeV, chi2/ndf=" << r.chi2_ndf << "\n";
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& ex) {
    std::cerr << "ERROR: " << ex.what() << "\n";
    return EXIT_FAILURE;
  }
}
