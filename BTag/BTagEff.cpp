#include "iostream"
#include "TClonesArray.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "TChain.h"
#include "classes/DelphesClasses.h"

using namespace std;

int coutBJet(TClonesArray *branchJet)
{
    int nBJet = 0;
    Jet *jet;

    for (int i = 0; i < branchJet->GetEntries(); i++)
    {
        jet = (Jet*) branchJet->At(i);
        if (jet->BTag == 1 && jet->PT > 20)
        {
            nBJet += 1;
        }
            
    }

    return nBJet;    
}

int main(int argc, char *argv[])
{
    // Usage: ./BTagEff /path/to/root/file

    TChain *chain = new TChain("Delphes");
    chain->Add(argv[1]);
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchJet = treeReader->UseBranch("Jet");

    int diBJet = 0, singleB = 0;

    for (int i = 0; i < treeReader->GetEntries(); i++)
    {
        treeReader->ReadEntry(i);
        int nBJet = coutBJet(branchJet);
        if (nBJet == 1)
        {
            singleB += 1;
        }
        else if (nBJet == 2)
        {
            diBJet += 1;
        }
    }
    cout << argv[1] << " " << "SingleB:" << singleB << " " << "diB:" << diBJet << " " << "BTagEff:" << (diBJet + 0.5 * singleB) / 10000 << endl;

    return 0;
}