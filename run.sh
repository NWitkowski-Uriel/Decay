#!/bin/bash
# run.sh – Build and run DeltaAnalysis executables
# Usage: ./run.sh [gui|compute|mt|rapidity|all]

set -e  # exit on error

# ------------------------------------------------------------
# 1. Check dependencies
# ------------------------------------------------------------
command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found"; exit 1; }
command -v root-config >/dev/null 2>&1 || { echo "ERROR: ROOT not found"; exit 1; }

# ------------------------------------------------------------
# 2. Set OpenMP environment (use all available cores)
# ------------------------------------------------------------
export OMP_NUM_THREADS=$(nproc)
export OMP_PROC_BIND=true
export OMP_PLACES=cores
echo "OpenMP: using $OMP_NUM_THREADS threads"

# ------------------------------------------------------------
# 3. Build the project
# ------------------------------------------------------------
BUILD_DIR="build"
if [ ! -d "$BUILD_DIR" ]; then
    mkdir "$BUILD_DIR"
fi
cd "$BUILD_DIR"
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# ------------------------------------------------------------
# 4. Copy experimental data (if present)
# ------------------------------------------------------------
if [ -d "../data" ]; then
    cp ../data/*.txt . 2>/dev/null || echo "No data files copied"
else
    echo "No data/ directory found – skipping data copy"
fi

# ------------------------------------------------------------
# 5. Determine what to run
# ------------------------------------------------------------
RUN_TARGET="gui"  # default
if [ $# -ge 1 ]; then
    RUN_TARGET="$1"
fi

case "$RUN_TARGET" in
    gui)
        echo "Starting GUI..."
        ./gui_analysis
        ;;
    compute)
        echo "Running compute_yields..."
        time ./compute_yields
        ;;
    mt)
        echo "Running plot_mt..."
        time ./plot_mt
        ;;
    rapidity)
        echo "Running plot_rapidity..."
        time ./plot_rapidity
        ;;
    all)
        echo "Running all calculation and plotting programs..."
        time ./compute_yields
        time ./plot_mt
        time ./plot_rapidity
        ;;
    *)
        echo "Usage: $0 [gui|compute|mt|rapidity|all]"
        echo "  gui       – start the Qt GUI (default)"
        echo "  compute   – run compute_yields only"
        echo "  mt        – run plot_mt only"
        echo "  rapidity  – run plot_rapidity only"
        echo "  all       – run compute_yields, then plot_mt, then plot_rapidity"
        exit 1
        ;;
esac

cd ..