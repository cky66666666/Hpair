cmake_minimum_required(VERSION 3.10)

project(eventAnalysis)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT CONFIG REQUIRED)

include_directories(/mnt/d/packages/ExRootAnalysis /mnt/d/packages/Delphes /mnt/d/work/Hpair/cpp/showerAnalysis/Tagger 
/mnt/d/work/Hpair/cpp/showerAnalysis)
add_subdirectory(Tagger)

add_executable(eventAnalysis eventAnalysis.cpp)

target_link_libraries(eventAnalysis PUBLIC ROOT::Core ROOT::Gpad ROOT::Graf3d ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore
ROOT::Matrix ROOT::MultiProc ROOT::Net ROOT::Physics ROOT::Postscript ROOT::RIO ROOT::ROOTDataFrame ROOT::ROOTVecOps
ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree /mnt/d/packages/Delphes/libDelphes.so /usr/local/lib/libfastjet.so 
/mnt/d/packages/ExRootAnalysis/libExRootAnalysis.so Tagger)