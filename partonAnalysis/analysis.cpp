#include "iostream"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "fastjet/ClusterSequence.hh"
using namespace std;

void jetAnalysis(TClonesArray *branchJet, ExRootTreeReader *treeReader){
    TH1I *histJet = new TH1I("jet_pt", "p_{T}^{j}",20, 0, 20);
    Jet *jet;
    TLorentzVector jetMomentum;
    vector<double> jetPT;
    for (int i = 0; i < treeReader->GetEntries(); i++)
    {
        treeReader->ReadEntry(i);
        /* for (int j = 0; j < branchJet->GetEntriesFast(); j++)
        {
            jet = (Jet*) branchJet->At(j);
            jetMomentum = jet->P4();
            jetPT.push_back(jetMomentum.Pt());
        } */
        histJet->Fill(branchJet->GetEntriesFast());
        /* jetPT.clear(); */
    }
    histJet->Draw();
}

void particleAnalysis(TClonesArray *branchParticle, ExRootTreeReader *treeReader){
    TH1I *histParticle = new TH1I("pdg", "pdg", 3, 20, 23);
    GenParticle *particle;
    for (int entry = 0; entry < treeReader->GetEntries(); entry++)
    {
        treeReader->ReadEntry(entry);
        for (int i = 0; i < branchParticle->GetEntriesFast(); i++)
        {
            particle = (GenParticle*) branchParticle->At(i);
            histParticle->Fill(particle->PID);
        }
    }
    histParticle->Draw();
}

int nBTagging(TClonesArray *branchJet){
    Jet *jet;
    int bTagging = 0;
    for (int entry = 0; entry < branchJet->GetEntriesFast(); entry++)
    {
        jet = (Jet*) branchJet->At(entry);
        if (jet->BTag == 1)
        {
            bTagging += 1;
        }
    }
    return bTagging;
}

void eventSelection(ExRootTreeReader *treeReader){
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
}


void analysis(){
    /* gROOT->ProcessLine(".include </mnt/d/work/Hpair/Delphes>");
    gROOT->ProcessLine(".include </mnt/d/work/Hpair/Delphes>/external"); */
    gSystem->Load("/mnt/d/work/Hpair/Delphes/libDelphes");
    TChain *chain = new TChain("Delphes");
    chain->Add("/mnt/d/work/Hpair/events/hhjHLLHC.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    Jet *jet;
    treeReader->ReadEntry(10);
    jet = (Jet*) branchJet->At(2);
    cout << "SqrtR:" << jet->PT << endl;
}

int main(){
    fastjet::ClusterSequence clus;
    analysis();
    return 0;
}