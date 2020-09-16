#include "iostream"

#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootConfReader.h"
using namespace std;

int main(int argc, char *argv[]){
    //Usage ./Tagger conf_file.tcl hepmc_file.hepmc 
    ExRootConfReader *confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);
    int maxEvent = confReader->GetInt("::MaxEvents", 0);
    int skipEvent = confReader->GetInt("::SkipEvents", 0);
    Delphes *modularDelphes = new Delphes("Delphes");
    return 0;
}
