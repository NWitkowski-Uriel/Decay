# Profiling & Benchmark Guide

## 1. Quick benchmark

```
./build/benchmark_performance
```

## 2. Wall-time measurement

```
/usr/bin/time -v ./compute_spectra
```

## 3. perf (CPU hotspots)

```
perf record -g ./compute_spectra
perf report
```

Look for:
- Integrate1D_high
- Integrate2D_high
- dN_dmt_* functions

## 4. Flamegraph (optional)

```
perf script | stackcollapse-perf.pl | flamegraph.pl > flame.svg
```

## 5. Scaling test

```
export OMP_NUM_THREADS=1
./benchmark_performance

export OMP_NUM_THREADS=8
./benchmark_performance
```

Compare scaling.

## Expected hotspots
- mt integrals
- BW/PS integration (if cache not used)

## Goal
After optimizations:
- BW/PS should disappear from hotspots
- Integrals dominate remaining cost
