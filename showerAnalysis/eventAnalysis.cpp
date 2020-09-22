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

#include "Tagger/BoostedHiggs.h"

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

vector<GenParticle*> getParton(TClonesArray *branchParticle)
{
    vector<GenParticle*> parton;
    GenParticle *particle;
    int nParticle = branchParticle->GetEntriesFast();
    for (int iParticle = 0; iParticle < nParticle; iParticle++)
    {
        particle = (GenParticle*) branchParticle->At(iParticle);
        if (abs(particle->PID) <= 5 || particle->PID == 21)
        {
            parton.push_back(particle);
        }
    }
    return parton;
}

int flavourAssociation(PseudoJet jet, vector<GenParticle*> parton)
{
    int flavour = -1;
    double deltaR = 0.5;
    TLorentzVector pJet;
    pJet.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.e());
    for (int i = 0; i < parton.size(); i++)
    {
        int pid = abs(parton[i]->PID);
        if (pid == 21) pid = 0;
        if (pJet.DeltaR(parton[i]->P4()) < deltaR && pid > flavour)
        {
            flavour = pid;
        }
    }
    return flavour;
}

vector<TLorentzVector> clusterRemmant(vector<PseudoJet> remmant, vector<GenParticle*> parton)
{
    ClusterSequence *sequence = new ClusterSequence(remmant, JetDefinition(antikt_algorithm, 0.4));
    vector<PseudoJet> remmantJet = sorted_by_pt(sequence->inclusive_jets(30));
    vector<PseudoJet> allJet = remmantJet;
    vector<int> higgsIndex = {};
    TLorentzVector softHiggs, hardJet;
    softHiggs.SetPxPyPzE(0, 0, 0, 0);
    hardJet.SetPxPyPzE(0, 0, 0 ,0);
    double deltaInvMass = 1000;

    int nJet = remmantJet.size();
    for (int i = 0; i < nJet; i++)
    {
        if (flavourAssociation(remmantJet[i], parton) != 5) continue;
        for (int j = i + 1; j < nJet; j++)
        {
            if (flavourAssociation(remmantJet[j], parton) != 5) continue;
            if (abs((remmantJet[i] + remmantJet[j]).m() - 125) < deltaInvMass)
            {
                deltaInvMass = abs((remmantJet[i] + remmantJet[j]).m() - 125);
                softHiggs.SetPxPyPzE((remmantJet[i] + remmantJet[j]).px(), (remmantJet[i] + remmantJet[j]).py(), (remmantJet[i] + remmantJet[j]).pz(), (remmantJet[i] + remmantJet[j]).e());
                higgsIndex = {i, j};
            }
            
        }
        
    }

    if (higgsIndex.size() == 2)
    {
        remmantJet.erase(remmantJet.begin() + higgsIndex[1]);
        remmantJet.erase(remmantJet.begin() + higgsIndex[0]);
    }

    if (remmantJet.size() > 0) hardJet.SetPxPyPzE(remmantJet[0].px(), remmantJet[0].py(), remmantJet[0].pz(), remmantJet[0].e());
    

    return {softHiggs, hardJet};
}

bool eventSelector(TLorentzVector hardHiggs, vector<TLorentzVector> remmantObject)
{
    bool status = true;
    if (hardHiggs.M() == 0 || remmantObject[0].M() == 0)
    {
        status = false;
    }
    if (hardHiggs.Pt() < 150 || abs(remmantObject[0].M() - 125) > 10 || remmantObject[1].Pt() < 200)
    {
        status = false;
    }
    return status;
}

int main(int argc, char *argv[])
{
    gSystem->Load("/mnt/d/work/Hpair/Delphes/libDelphes");
    
    TChain *chain = new TChain("Delphes");
    chain->Add("/mnt/d/work/Hpair/events/root/bkg4b.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");

    TLorentzVector hardHiggs;
    vector<TLorentzVector> remmantObject;

    vector<PseudoJet> finalState;
    vector<GenParticle*> parton;

    int nEvent = treeReader->GetEntries();
    int n = 0, m = 0;
    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        if(trigger(branchJet))
        {
            finalState = getFinalState(branchTower);
            parton = getParton(branchParticle);
            BoostedHiggs *boostedHiggs = new BoostedHiggs(finalState, parton, 150, 110, 1.5, 0.3);
            boostedHiggs->process();
            if ((boostedHiggs->boostedHiggs.size()) < 2) continue;
            PseudoJet tmp = (boostedHiggs->boostedHiggs)[0] + (boostedHiggs->boostedHiggs)[1];
            hardHiggs.SetPxPyPzE(tmp.px(), tmp.py(), tmp.pz(), tmp.e());
            remmantObject = clusterRemmant(boostedHiggs->remmant, parton);
            if (eventSelector(hardHiggs, remmantObject))
            {
                m += 1;
            }
            delete boostedHiggs;
        }
    }
    cout << "mHiggs" << m << endl;
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