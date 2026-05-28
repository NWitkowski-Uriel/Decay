cmake_minimum_required(VERSION 3.16)
project(DeltaRoot3D VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree Graf Gpad MathCore)
find_package(OpenMP)

add_library(delta3d STATIC
  src/Physics.cpp
)

target_include_directories(delta3d PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${ROOT_INCLUDE_DIRS}
)

target_compile_options(delta3d PRIVATE -O3 -Wall -Wextra -Wpedantic)
if(OpenMP_CXX_FOUND)
  target_link_libraries(delta3d PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(delta3d PUBLIC DELTA3D_USE_OPENMP=1)
endif()

target_link_libraries(delta3d PUBLIC ${ROOT_LIBRARIES})

add_executable(delta3d_run apps/delta3d_run.cpp)
target_link_libraries(delta3d_run PRIVATE delta3d)

install(TARGETS delta3d_run RUNTIME DESTINATION bin)
