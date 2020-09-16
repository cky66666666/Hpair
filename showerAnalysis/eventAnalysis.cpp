#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TObject.h"

#include "iostream"
#include "vector"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

bool trigger(TClonesArray *branchJet)
{
    Jet *jet;
    vector<double> jetPtVec; 
    int nJet = branchJet->GetEntriesFast();
    for (int iJet = 0; iJet < nJet; iJet++)
    {
        jet = (Jet*) branchJet->At(iJet);
        if (!(jet->IsTotalParticle))
        {
            jetPtVec.push_back(jet->PT);
        }
    }
    if (jetPtVec.size() < 4)
    {
        return false;
    }
    else
    {
        sort(jetPtVec.begin(), jetPtVec.end());
        reverse(jetPtVec.begin(), jetPtVec.end());
        if (jetPtVec[0] > 120 && jetPtVec[1] > 100 && jetPtVec[2] > 70 && jetPtVec[3] > 40)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

vector<PseudoJet> getFinalState(TClonesArray *branchTower)
{
    vector<PseudoJet> finalState;
    Tower *tower;
    GenParticle *particle;
    TLorentzVector momentum;
    momentum.SetPxPyPzE(0, 0, 0, 0);
    PseudoJet tmp;
    int nTower = branchTower->GetEntriesFast();
    int userIndex = 0;
    for (int iTower = 0; iTower < nTower; iTower++)
    {
        tower = (Tower*) branchTower->At(iTower);
        int nParticle = tower->Particles.GetEntriesFast();
        for (int iParticle = 0; iParticle < nParticle; iParticle++)
        {
            if (tower->Particles.At(iParticle)->IsA() == GenParticle::Class())
            {
                particle = static_cast<GenParticle *> (tower->Particles.At(iParticle));
                momentum = particle->P4();
                tmp = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
                tmp.set_user_index(userIndex);
                finalState.push_back(tmp);
                userIndex += 1; 
            }
        }
    }
    
    return finalState;
}

int main(int argc, char *argv[])
{
    gSystem->Load("/mnt/d/work/Hpair/Delphes/libDelphes");
    
    TChain *chain = new TChain("Delphes");
    chain->Add("/mnt/d/work/Hpair/events/root/hhj1.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");

    vector<PseudoJet> finalState;
    int nEvent = treeReader->GetEntries();
    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        finalState = getFinalState(branchTower);
        int n = 0;
        for (int i = 0; i < branchTower->GetEntriesFast(); i++)
        {
            Tower *tower = (Tower*) branchTower->At(i);
            n += tower->Particles.GetEntriesFast();
        }
        
        if ((n - finalState.size()) != 0) cout << finalState.size() << endl;
    }
    /* TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    treeReader->ReadEntry(1);
    Jet *jet;
    GenParticle *particle;
    int n = 0;
    for (int i = 0; i < branchJet->GetEntries(); i++)
    {
        jet = (Jet*) branchJet->At(i);
        cout << jet->IsTotalParticle << endl;
        if (jet->IsTotalParticle)
        {
            cout << (jet->Constituents).GetEntries() << endl;
        }
        else
        {
            n += (jet->Constituents).GetEntries();
        }
    }
    cout << n << endl; */
    return 0;
}