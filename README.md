# DeltaAnalysis

Repozytorium zawiera dawne wersje programów CERN ROOT oraz nowy, czysty moduł `root3d/`, przygotowany jako baza pod aktualne notatniki Wolframa w wersji 3D.

## Zalecany nowy moduł: `root3d/`

Nowa implementacja nie nadpisuje starych plików.  Ma osobny system budowania i osobny katalog wyników.

```bash
cmake -S root3d -B build-root3d
cmake --build build-root3d -j
./build-root3d/delta3d_run --fast
```

Pełniejszy opis znajduje się w:

```text
root3d/README.md
```

## Co zawiera `root3d/`

- krotności pierwotne, poprawki i sumy,
- modele masy: Dirac, Breit-Wigner, Phase Shift,
- kanały `Delta++`, `Delta+`, `Delta0`, `Delta-` z czynnikami `1`, `2/3`, `1/3`,
- testy normalizacji gałęzi w stylu aktualnego notebooka Yield Tests,
- widma `mT` i rapidity,
- eksport CSV, ROOT `TTree`/`TGraph` i wykresy PNG,
- bezpieczniejszą parallelizację po zewnętrznych punktach siatki, nie przez rozbijanie pojedynczego całkowania.

## Status naukowy

Warstwa krotności, kanałów, modeli masy i eksportu jest przygotowana pod aktualny workflow.  Jedyny element oznaczony jako wymagający bezpośredniego portu z finalnego notebooka to dokładny kształt widma poprawki córki z pełnym całkowaniem q-k-r.  Obecnie kod zachowuje poprawną normalizację poprawki i daje kontrolny kształt termiczny, ale finalne publikacyjne porównania wymagają wklejenia dokładnego kernela z Wolframa do funkcji wskazanych w `root3d/src/Physics.cpp`.

## Starszy kod

Starsze pliki ROOT/GUI pozostają w historii repozytorium.  Nowe prace najlepiej prowadzić w `root3d/`, aby uniknąć mieszania nieaktualnego pipeline'u z aktualną analizą 3D.
