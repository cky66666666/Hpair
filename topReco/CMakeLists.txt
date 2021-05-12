cmake_minimum_required(VERSION 3.10)

project(fit)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT CONFIG REQUIRED)

include_directories(/home/E/chaikangyu/packages/include /home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/Tagger
/home/E/chaikangyu/packages/MG5_aMC_v2_6_7/Delphes /home/E/chaikangyu/work/Hpair/cpp/showerAnalysis /home/E/chaikangyu/packages/ExRootAnalysis)

add_executable(fit fit.cpp)

target_link_libraries(fit PUBLIC ROOT::Core ROOT::Gpad ROOT::Graf3d ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore
ROOT::Matrix ROOT::MultiProc ROOT::Net ROOT::Physics ROOT::Postscript ROOT::RIO ROOT::ROOTDataFrame ROOT::ROOTVecOps
ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree /home/E/chaikangyu/packages/MG5_aMC_v2_6_7/Delphes/libDelphes.so /home/E/chaikangyu/packages/lib/libfastjet.so 
/home/E/chaikangyu/packages/ExRootAnalysis/libExRootAnalysis.so PUBLIC)