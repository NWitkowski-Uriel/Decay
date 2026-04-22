# CMakeLists.txt for Delta Analysis Project
# Requires CMake 3.10 or higher
cmake_minimum_required(VERSION 3.10)

project(DeltaAnalysis VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable automatic moc and uic for Qt
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree Graf Gpad)
find_package(OpenMP REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Widgets)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
    ${Qt5Widgets_INCLUDE_DIRS}
)

# Static library with core functionality
set(LIB_SOURCES
    src/MathUtils.cpp
    src/Thermodynamics.cpp
    src/DecayFunctions.cpp
    src/Plotting.cpp
    src/RadialGrid.cpp
)

add_library(DeltaLib STATIC ${LIB_SOURCES})
target_link_libraries(DeltaLib ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_options(DeltaLib PRIVATE -O3 -march=native -ffast-math)

# ----------------------------------------------------------------------
# Executables
# ----------------------------------------------------------------------

# Original compute_yields (kept for compatibility, but not used by GUI)
add_executable(compute_yields apps/compute_yields.cpp)
target_link_libraries(compute_yields DeltaLib)
target_compile_options(compute_yields PRIVATE -O3 -march=native -ffast-math)

# New specialised executables
add_executable(compute_primordial apps/compute_primordial.cpp)
target_link_libraries(compute_primordial DeltaLib)
target_compile_options(compute_primordial PRIVATE -O3 -march=native -ffast-math)

add_executable(compute_delta_decays apps/compute_delta_decays.cpp)
target_link_libraries(compute_delta_decays DeltaLib)
target_compile_options(compute_delta_decays PRIVATE -O3 -march=native -ffast-math)

add_executable(compute_total apps/compute_total.cpp)
target_link_libraries(compute_total DeltaLib)
target_compile_options(compute_total PRIVATE -O3 -march=native -ffast-math)

# Plotting programs
add_executable(plot_mt apps/plot_mt_distributions.cpp)
target_link_libraries(plot_mt DeltaLib)
target_compile_options(plot_mt PRIVATE -O3 -march=native -ffast-math)

add_executable(plot_rapidity apps/plot_rapidity.cpp)
target_link_libraries(plot_rapidity DeltaLib)
target_compile_options(plot_rapidity PRIVATE -O3 -march=native -ffast-math)

# GUI
add_executable(gui_analysis
    gui/main.cpp
    gui/mainwindow.cpp
)
target_link_libraries(gui_analysis
    DeltaLib
    ${ROOT_LIBRARIES}
    OpenMP::OpenMP_CXX
    Qt5::Widgets
)
target_compile_options(gui_analysis PRIVATE -O3 -march=native -ffast-math)

# ----------------------------------------------------------------------
# Messages
# ----------------------------------------------------------------------
message(STATUS "DeltaAnalysis project configured.")
message(STATUS "ROOT version: ${ROOT_VERSION}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "Qt5 version: ${Qt5_VERSION}")