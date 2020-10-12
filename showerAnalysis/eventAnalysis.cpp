#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "THStack.h"
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

struct Remmant
{
    Remmant()
    {
        hardJet.SetPxPyPzE(0, 0, 0, 0);
        softHiggs.SetPxPyPzE(0, 0, 0, 0);
    }
    TLorentzVector hardJet;
    TLorentzVector softHiggs;
};

struct myEvent
{
    myEvent()
    {
        hardJet.SetPxPyPzE(0, 0, 0, 0);
        hardHiggs.SetPxPyPzE(0, 0, 0, 0);
        softHiggs.SetPxPyPzE(0, 0, 0, 0);
    }
    TLorentzVector hardJet;
    TLorentzVector hardHiggs;
    TLorentzVector softHiggs;
};



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

Remmant clusterRemmant(vector<PseudoJet> remmant, vector<GenParticle*> parton, int type)
{
    //type=1: no higgs tagger
    //type=2: higgs tagger for h->2b
    //type=3: higgs tagger for h->2a
    ClusterSequence *sequence = new ClusterSequence(remmant, JetDefinition(antikt_algorithm, 0.4));
    vector<PseudoJet> remmantJet = sorted_by_pt(sequence->inclusive_jets(30));
    
    Remmant remmantObject;

    if (type == 1)
    {
        TLorentzVector hardJet;
        if (remmantJet.size() > 0) hardJet.SetPxPyPzE(remmantJet[0].px(), remmantJet[0].py(), remmantJet[0].pz(), remmantJet[0].e());
        remmantObject.hardJet = hardJet;
        delete sequence;
        return remmantObject;
    }
    else if (type == 2)
    {
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
        
        remmantObject.hardJet = hardJet;
        remmantObject.softHiggs = softHiggs;

        delete sequence;
        return remmantObject;
    }
    else
    {
        
        delete sequence;
        return remmantObject;
    }
    
}

bool eventSelector(myEvent event)
{
    bool status = true;
    if (event.hardHiggs.M() == 0)
    {
        status = false;
    }
    if (abs(event.softHiggs.M() - 125) > 20 || event.hardJet.Pt() < 200)
    {
        status = false;
    }
    return status;
}

double analyse4b(TClonesArray *branchJet, TClonesArray *branchParticle, TClonesArray *branchTower, int type)
{
    //type=1: single boosted higgs
    //type=2: double boosted higgs
    //type=3: collinear higgs pair

    vector<PseudoJet> finalState;
    vector<GenParticle*> parton;
    BoostedHiggs *boostedHiggs = new BoostedHiggs();

    PseudoJet tmp;
    Remmant remmantObject;
    myEvent event;
    
    if(trigger(branchJet))
    {
        finalState = getFinalState(branchTower);
        parton = getParton(branchParticle);

        boostedHiggs->init(finalState, parton, 150, 110, 1.5, 0.3);
        if (type == 1)
        {
            boostedHiggs->process(1);
            if (boostedHiggs->boostedHiggs.size() < 1)
            {
                boostedHiggs->clear();
                delete boostedHiggs;
                return 0;
            }
            tmp = (boostedHiggs->boostedHiggs)[0];
            event.hardHiggs.SetPxPyPzE(tmp.px(), tmp.py(), tmp.pz(), tmp.e());

            remmantObject = clusterRemmant(boostedHiggs->remmant, parton, 2);
            event.hardJet = remmantObject.hardJet;
            event.softHiggs = remmantObject.softHiggs;
        }
        else if (type == 2 || type == 3)
        {
            boostedHiggs->process(type);
            if (boostedHiggs->boostedHiggs.size() < 2)
            {
                boostedHiggs->clear();
                delete boostedHiggs;
                return 0;
            }
            tmp = (boostedHiggs->boostedHiggs)[0];
            event.hardHiggs.SetPxPyPzE(tmp.px(), tmp.py(), tmp.pz(), tmp.e());
            tmp = (boostedHiggs->boostedHiggs)[1];
            event.softHiggs.SetPxPyPzE(tmp.px(), tmp.py(), tmp.pz(), tmp.e());

            remmantObject = clusterRemmant(boostedHiggs->remmant, parton, 1);
            event.hardJet = remmantObject.hardJet;
        }
        if (eventSelector(event))
        {
            boostedHiggs->clear();
            delete boostedHiggs;
            return (event.hardHiggs + event.softHiggs).M();
            //cout << (event.hardHiggs + event.softHiggs).M() << endl;
        }
        else
        {
            delete boostedHiggs;
            return 0;
        }
    }
    else
    {
        delete boostedHiggs;
        return 0;
    }
}

int main(int argc, char *argv[])
{
    // Usage: ./eventAnalysis /path/to/hist.root /path/to/inputfile histogram name

    gSystem->Load("/mnt/d/work/Hpair/Delphes/libDelphes");
    
    TH1D *hist = new TH1D(argv[3], argv[3], 50, 250, 1000);
    TFile *f = new TFile(argv[1], "RECREATE");

    TChain *chain = new TChain("Delphes");
    chain->Add(argv[2]);
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");

    //cout << treeReader->GetEntries() << endl;
    int nEvent = treeReader->GetEntries();
    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        //cout << iEvent << endl;
        double inv = analyse4b(branchJet, branchParticle, branchTower, 1);
        if (inv != 0)
        {
            hist->Fill(inv);
        }
        treeReader->Clear();
    }
    
    hist->Write();    
    f->Close();

    return 0;
}