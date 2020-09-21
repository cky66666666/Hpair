#include "iostream"
#include "vector"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootClasses.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "algorithm"

using namespace std;

vector<TRootLHEFParticle*> findParticles(TClonesArray *branchParticles){
    vector<TRootLHEFParticle*> particleList;
    TRootLHEFParticle *particle;
    for (int i = 0; i < branchParticles->GetEntries(); i++)
    {
        particle = (TRootLHEFParticle*) branchParticles->At(i);
        if (particle->Status == 1)
        {
            particleList.push_back(particle);
        }
    }
    return particleList;
}

int maxPt(vector<TLorentzVector> object){
    vector<double> pt;
    double maxValue;
    int maxPoint;
    for (int i = 0; i < object.size(); i++)
    {
        pt.push_back(object[i].Pt());
    }
    maxValue = pt[0];
    maxPoint = 0;
    for (int j = 0; j < pt.size(); j++)
    {
        if (pt[j] > maxValue)
        {
            maxValue = pt[j];
            maxPoint = j;
        }
    }
    return maxPoint;
}

int minPt(vector<TLorentzVector> object){
    vector<double> pt;
    double minValue;
    int minPoint;
    for (int i = 0; i < object.size(); i++)
    {
        pt.push_back(object[i].Pt());
    }
    minValue = pt[0];
    minPoint = 0;
    for (int j = 0; j < pt.size(); j++)
    {
        if (pt[j] < minValue)
        {
            minValue = pt[j];
            minPoint = j;
        }
    }
    return minPoint;
}

TH1D* drawHist(ExRootTreeReader *treeReader, int legend){
    TClonesArray *branchParticles = treeReader->UseBranch("Particle");
    TRootLHEFParticle *particle;
    vector<TLorentzVector> higgs;
    TLorentzVector pHiggsPair, tmp;
    int maxPtPoint, minPtPoint;
    char histName[10];
    sprintf(histName, "deltaLam=%d", legend - 6);
    TH1D *histInvMass = new TH1D(histName, histName, 50, 250, 1000);
    for (int iEvent = 0; iEvent < treeReader->GetEntries(); iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        pHiggsPair.SetPxPyPzE(0, 0, 0, 0);
        for (int iParticle = 0; iParticle < branchParticles->GetEntries(); iParticle++)
        {
            particle = (TRootLHEFParticle*) branchParticles->At(iParticle);
            if (particle->PID == 25)
            {
                tmp.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
                higgs.push_back(tmp);
            }
        }
        //int maxPoint = maxPt(higgs);
        histInvMass->Fill((higgs[0] + higgs[1]).M());
        higgs.clear();
    }
    return histInvMass;
}

TLorentzVector hardHiggsMomentum(TClonesArray *branchParticle){
    TLorentzVector pHard;
    TRootLHEFParticle *particle;
    pHard.SetPxPyPzE(0, 0, 0, 0);
    for (int i = 0; i < branchParticle->GetEntries(); i++)
    {
        particle = (TRootLHEFParticle*) branchParticle->At(i);
        if (particle->PID == 25 && (particle->E) > pHard.E())
        {
            pHard.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
        }
    }
    return pHard;
}

void deltaRFilter(ExRootTreeReader *treeReader){
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TLorentzVector pHiggs, partonMomentum;
    vector<TRootLHEFParticle*> candidate;
    int n = 0;
    for (int iEvent = 0; iEvent < treeReader->GetEntries(); iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        pHiggs = hardHiggsMomentum(branchParticle);
        for (int iParticle = 0; iParticle < branchParticle->GetEntries(); iParticle++)
        {
            TRootLHEFParticle *particle = (TRootLHEFParticle*) branchParticle->At(iParticle);
            if (particle->Status == 1)
            {
                partonMomentum.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
                if (pHiggs.DeltaR(partonMomentum) < 1.5)
                {
                    candidate.push_back(particle);
                }
            }
        }
        for (int i = 0; i < candidate.size(); i++)
        {
            if (abs(candidate[i]->PID) != 5)
            {
                n += 1;
            }
            
        }
        
        candidate.clear();
    }
    cout << n << endl;
}

bool eventSelector(TRootLHEFParticle *particle)
{
    
}

int main(int argc, char *argv[]){

    char treeName[] = "bkg4b";
    TChain *chain = new TChain(treeName);
    chain->Add("../bkg4b.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

    TClonesArray *branchParticle;
    branchParticle = treeReader->UseBranch("Particle");

    TRootLHEFParticle *particle;

    int nEvent = treeReader->GetEntries();
    int n = 0;
    double minPt;
    cout << nEvent << endl;
    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        for (int iParticle = 0; iParticle < branchParticle->GetEntriesFast(); iParticle++)
        {
            particle = (TRootLHEFParticle*) branchParticle->At(iParticle);
            if ((particle->Status == 1) /* && (abs(particle->PID) != 5) */ && (particle->PT > 100))
            {
                n += 1;
            }
            
        }
        
    }
    
    cout << n << endl;

    /* vector<int> nTree = {6, 7, 8, 9, 10};
    THStack *stackInvMass = new THStack("InvMass", "InvMass");
    for (int i = 0; i < nTree.size(); i++)
    {
        char treeName[10];
        sprintf(treeName, "LHEF%d", nTree[i]);
        TChain *chain = new TChain(treeName);
        chain->Add(argv[1]);
        ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
        stackInvMass->Add(drawHist(treeReader, nTree[i]));
    } */

    /* TFile f("../hist.root", "RECREATE");
    stackInvMass->Write();
    f.Close(); */
    return 0;
}