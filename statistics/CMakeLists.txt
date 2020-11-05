cmake_minimum_required(VERSION 3.10)

project(statistics)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT CONFIG REQUIRED)

include_directories(/mnt/d/work/Hpair/cpp/statistics/include /mnt/d/work/Hpair/cpp/statistics)
add_subdirectory(include)

add_executable(statistics statistics.cpp)

target_link_libraries(statistics PUBLIC ROOT::Core ROOT::Gpad ROOT::Graf3d ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore
ROOT::Matrix ROOT::MultiProc ROOT::Net ROOT::Physics ROOT::Postscript ROOT::RIO ROOT::ROOTDataFrame ROOT::ROOTVecOps
ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree include PUBLIC)