cmake_minimum_required(VERSION 3.10)

project(partonAnalysis)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT CONFIG REQUIRED)

include_directories(/home/E/chaikangyu/packages/MG5_aMC_v2_6_7/ExRootAnalysis /home/E/chaikangyu/packages/MG5_aMC_v2_6_7/Delphes)

add_executable(partonAnalysis partonAnalysis.cpp)

target_link_libraries(partonAnalysis PUBLIC ROOT::Core ROOT::Gpad ROOT::Graf3d ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore
ROOT::Matrix ROOT::MultiProc ROOT::Net ROOT::Physics ROOT::Postscript ROOT::RIO ROOT::ROOTDataFrame ROOT::ROOTVecOps
ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree /home/E/chaikangyu/packages/MG5_aMC_v2_6_7/ExRootAnalysis/libExRootAnalysis.so)