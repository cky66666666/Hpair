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
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

    int diBJet = 0, singleB = 0;
    int diPhoton = 0, singlePhoton = 0;

    for (int i = 0; i < treeReader->GetEntries(); i++)
    {
        treeReader->ReadEntry(i);
        int nBJet = coutBJet(branchJet);
        int nPhoton = branchPhoton->GetEntries();
        if (nPhoton == 1)
        {
            singlePhoton += 1;
        }
        else if (nPhoton == 2)
        {
            diPhoton += 1;
        }
    }
    cout << argv[1] << " " << "SinglePhoton:" << singlePhoton << " " << "diPhoton:" << diPhoton << " " << "PhotonEff:" << (diPhoton + 0.5 * singlePhoton) / 10000 << endl;

    return 0;
}