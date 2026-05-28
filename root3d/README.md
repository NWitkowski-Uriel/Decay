# Delta ROOT 3D rewrite

This directory is a clean CERN ROOT/C++ rewrite of the current 3D Wolfram workflow for Delta resonance decay studies.

It is intentionally placed in `root3d/` instead of replacing the old project directly.  The older ROOT code and repository history stay available, while this module gives a reproducible baseline that can be extended notebook-by-notebook.

## Scope

Implemented now:

- 3D thermal/blast-wave style primordial kernel,
- stable particles: `p`, `n`, `pi+`, `pi0`, `pi-`,
- Delta parents: `Delta++`, `Delta+`, `Delta0`, `Delta-`,
- decay channel bookkeeping with branch factors `1`, `2/3`, `1/3`,
- mass models: `Dirac`, relativistic `BW`, and Pok-Man-Lo-style `PS`,
- normalized mass distributions on `[mN + mpi, 1500 MeV]`,
- multiplicity tables: primordial, correction, total,
- branch normalization tests analogous to Wolfram Yield Tests,
- `mT` and rapidity spectra for each stable particle, parameter set, and mass model,
- ROOT output with `TTree` and `TGraph` objects,
- PNG plot export.

Important current limitation:

- The daughter correction spectrum shape is normalized correctly to the decay correction yield, but its shape currently uses a controlled thermal proxy.  The exact Wolfram q-k-r daughter convolution must be copied into `correction_shape_mt` / `correction_shape_y` in `src/Physics.cpp` once the final current notebook cell is available.  The bookkeeping and normalization layer are already prepared for that replacement.

## Build

From the repository root:

```bash
cmake -S root3d -B build-root3d
cmake --build build-root3d -j
```

If ROOT is not on `PATH`, source your ROOT installation first, for example:

```bash
source /path/to/root/bin/thisroot.sh
```

For the user's known FairSoft-style installation this may be similar to:

```bash
source ~/fairsoft_nov22p4_root6/installation/bin/thisroot.sh
```

## Run

Smoke test:

```bash
./build-root3d/delta3d_run --fast
```

Default run:

```bash
./build-root3d/delta3d_run
```

Higher-resolution run:

```bash
./build-root3d/delta3d_run --production
```

Serial run, useful when validating numerical reproducibility:

```bash
./build-root3d/delta3d_run --serial
```

Custom output directory:

```bash
./build-root3d/delta3d_run --output output/my_run
```

## Outputs

Default output directory: `output/root3d/`.

Files:

- `multiplicities.csv` — columns: `parameter_set,mass_model,particle,primordial,correction,total`,
- `yield_tests.csv` — branch tests in the style of the current Yield Tests notebook,
- `delta3d.root` — ROOT file containing:
  - tree `multiplicities`,
  - directory `spectra/<parameter>_<model>_<particle>/`,
  - graphs `mt_primordial`, `mt_correction`, `mt_total`, `y_primordial`, `y_correction`, `y_total`,
- `plots/*.png` — quick-look plots.

## Where to paste fitted notebook values

Edit `DefaultRunConfig()` in `src/Physics.cpp`:

```cpp
ModelConfig fitted = vacuum;
fitted.tag = "fitted";
fitted.m_delta = 1232.0;
fitted.gamma_delta = 117.0;
```

Edit phase-shift parameters in `include/Delta3D/Config.h` or assign them per parameter set in `DefaultRunConfig()`:

```cpp
double alpha0_ps = ...;
double c1_ps = ...;
double c2_ps = ...;
```

These are deliberately explicit because the final values should be synchronized with the Wolfram PairFit and Distributions notebooks.

## Numerical strategy

The old attempt parallelized too deeply and could split single expensive spectra across kernels in ways that made debugging difficult.  This rewrite parallelizes outer grid loops only:

- each `mT` point is independent,
- each rapidity point is independent,
- low-level integration routines use fixed-grid Simpson/midpoint quadrature for deterministic behavior.

This mirrors the safer strategy discussed for Wolfram: one spectrum/table point is a unit of parallel work; the inner physics integral remains local and deterministic.

## Next porting step

The next scientifically important step is to replace the current proxy correction shape with the exact Wolfram daughter integrals:

- Dirac daughter kernel with q-k-r integration,
- BW/PS daughter kernel with outer M integration,
- same limits `km`, `kp`, `PStar`, and branch factors as in the notebooks.

The functions already exposed for this are:

```cpp
double PStar(double M, double m1, double m2);
double KMinus(double q, double M, double m1, double m2);
double KPlus(double q, double M, double m1, double m2);
double SpectralWeight(MassModel model, double M, const ModelConfig& cfg, double m1, double m2);
```
