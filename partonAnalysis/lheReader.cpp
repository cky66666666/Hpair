#include <stdexcept>
#include <iostream>
#include <sstream>

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootLHEFReader.h"
#include "TFile.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

int main(int argc, char *argv[]){
    //Usage: ./lheReader /path/to/output.root treename /path/to/input.lhe
    
    FILE *inputfile;
    TFile *outputfile;
    outputfile = TFile::Open(argv[1], "RECREATE");

    ExRootLHEFReader *reader = new ExRootLHEFReader;
    ExRootTreeWriter *treeWriter = new ExRootTreeWriter(outputfile, argv[2]);
    ExRootTreeBranch *branchEvent = 0, *branchRwgt = 0, *branchParticle = 0;
    branchEvent = treeWriter->NewBranch("Event", TRootLHEFEvent::Class());
    branchRwgt = treeWriter->NewBranch("Rwgt", TRootWeight::Class());
    branchParticle = treeWriter->NewBranch("Particle", TRootLHEFParticle::Class());
    
    inputfile = fopen(argv[3], "r");
    fseek(inputfile, 0L, SEEK_END);
    int length = ftello(inputfile);
    fseek(inputfile, 0L, SEEK_SET);

    reader->SetInputFile(inputfile);
    int eventCounter = 0;
    treeWriter->Clear();
    reader->Clear();
    while(reader->ReadBlock(branchParticle))
    {
        if(reader->EventReady())
        {
        ++eventCounter;

        reader->AnalyzeEvent(branchEvent, eventCounter);
        reader->AnalyzeRwgt(branchRwgt);

        treeWriter->Fill();

        treeWriter->Clear();

        reader->Clear();
        }
    }
    fseek(inputfile, 0L, SEEK_END);
    fclose(inputfile);
    treeWriter->Write();
    delete reader;
    delete treeWriter;

    return 0;
}