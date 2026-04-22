# CMakeLists.txt for Delta Analysis Project
cmake_minimum_required(VERSION 3.10)

project(DeltaAnalysis VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# New architecture apps
add_executable(compute_spectra apps/compute_spectra.cpp)
target_link_libraries(compute_spectra DeltaLib)

add_executable(compute_yields_refactored apps/compute_yields_refactored.cpp)
target_link_libraries(compute_yields_refactored DeltaLib)

add_executable(compare_models apps/compare_models.cpp)
target_link_libraries(compare_models DeltaLib)

add_executable(export_publication_plots apps/export_publication_plots.cpp)
target_link_libraries(export_publication_plots DeltaLib)

add_executable(test_decay_multiplicities apps/test_decay_multiplicities.cpp)
target_link_libraries(test_decay_multiplicities DeltaLib)

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

message(STATUS "DeltaAnalysis project configured.")
