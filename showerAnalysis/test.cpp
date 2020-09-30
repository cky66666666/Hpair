#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "iostream"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "TClonesArray.h"

using namespace std;

void test()
{
    TFile f("../../events/root/hhj1.root");
    TTree *t;
    TBranch *branch;
    int len;
    float pt[100];
    gSystem->Load("/mnt/d/packages/Delphes/libDelphes");
    t = (TTree*) f.Get("Delphes");

    branch = t->GetBranch("Particle.PT");
    t->SetMakeClass(1);
    t->SetBranchAddress("Particle", &len);

    branch->SetAddress(&pt);
    int n = t->GetEntries();
    cout << n << endl;
    for (int i = 0; i < n; i++)
    {
        t->GetEntry(i);
        //cout << branch->GetEntries() << endl;
    }
    
}