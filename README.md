# DeltaAnalysis (CERN ROOT, bez GUI)

Repozytorium zawiera implementację analizy z notatników Wolfram (`step3D` stabilny i `step4C` parallel) w postaci programów C++/ROOT.

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Główny pipeline (obliczenia + tabele + wykresy)

```bash
./build/compute_spectra
```

Wyniki:
- `output/step4c_analysis.root` – drzewo `multiplicity` z tabelą wyników,
- `output/multiplicity_table.csv` – tabela CSV,
- `output/step4c_*_mt.pdf/png` – wykresy mT dla wszystkich stabilnych cząstek.

## Opcjonalnie GUI

GUI jest domyślnie wyłączone. Aby je budować:

```bash
cmake -S . -B build -DBUILD_GUI=ON
cmake --build build -j
```
