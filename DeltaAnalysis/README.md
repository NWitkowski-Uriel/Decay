# DeltaAnalysis (ROOT + Qt GUI)

Aplikacja odtwarza fizykę z notatników Wolfram i umożliwia uruchamianie obliczeń z poziomu GUI:

1. **widma pierwotne** (`compute_primordial`),
2. **wkłady z rozpadów Δ** (`compute_delta_decays`, modele Dirac/Breit–Wigner/Phase Shift),
3. **sumy i normalizacja** (`compute_total`),
4. **wykresy mT i y** (`plot_mt`, `plot_rapidity`).

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Uruchomienie GUI

```bash
cd build
./gui_analysis
```

GUI uruchamia etapy obliczeń i zapisuje wyniki do plików ROOT (`primordial.root`, `delta.root`, `total.root`) oraz wykresy do katalogu `output/`.

## Uruchomienia CLI

### 1) Widma pierwotne

```bash
./compute_primordial
```

### 2) Rozpady Δ

```bash
./compute_delta_decays --particles=proton,piplus --model=dirac
# --model: dirac | bw | ps | all
```

### 3) Suma + skala

```bash
./compute_total
```

### 4) Wykresy

```bash
./plot_mt --distributions=mt,ratio --particles=proton,piplus
./plot_rapidity --distributions=rapidity,ratio --particles=proton,piplus
```

## Uwaga o zgodności z notatnikami

Kod `compute_delta_decays` korzysta z pełnego zestawu funkcji dla wszystkich cząstek (`p, n, pi+, pi-, pi0`) oraz trzech modeli masy rezonansu, zgodnie z implementacją w bibliotece `DecayFunctions`.
