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

double deltaR(PseudoJet j1, PseudoJet j2)
{
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(j1.px(), j1.py(), j1.pz(), j1.e());
    p2.SetPxPyPzE(j2.px(), j2.py(), j2.pz(), j2.e());
    return p1.DeltaR(p2);
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

vector<PseudoJet> findHiggs(vector<vector<PseudoJet>> higgsCandidate, vector<GenParticle*> parton)
{
    double deltaInvMass, mH = 125.0, tmp;
    vector<PseudoJet> candidate = {}, higgs;
    TLorentzVector p1, p2;
    //cout << "higgs" << " " << higgsCandidate.size() << endl;
    deltaInvMass = 1000;
    for (int i = 0; i < higgsCandidate.size(); i++)
    {
        tmp = 1000;
        for (int j = 0; j < higgsCandidate[i].size(); j++)
        {
            for (int k = j + 1; k < higgsCandidate[i].size(); k++)
            {
                p1.SetPxPyPzE(higgsCandidate[i][j].px(), higgsCandidate[i][j].py(), higgsCandidate[i][j].pz(), higgsCandidate[i][j].e());
                p2.SetPxPyPzE(higgsCandidate[i][k].px(), higgsCandidate[i][k].py(), higgsCandidate[i][k].pz(), higgsCandidate[i][k].e());
                p1 = p1 + p2;
                if (abs(p1.M() - mH) < tmp )
                {
                    tmp = abs(p1.M() - mH);
                    candidate = {higgsCandidate[i][j], higgsCandidate[i][k]}; 
                }
                
            }
            
        }
        if (candidate.size() < 2) continue;
        if (flavourAssociation(candidate[0], parton) == 5 && flavourAssociation(candidate[1], parton) == 5)
        {
            if (tmp < deltaInvMass)
            {
                deltaInvMass = tmp;
                higgs = candidate;
            }
            
        }
        
    }
    return higgs;
}

vector<PseudoJet> higgsTagger(vector<PseudoJet> finalState, vector<GenParticle*> parton)
{
    ClusterSequence *sequence = new ClusterSequence(finalState, JetDefinition(cambridge_algorithm, 1.5));
    vector<PseudoJet> fatjet = sorted_by_pt(sequence->inclusive_jets(150.0));
    vector<PseudoJet> mother, subStructure, tmp;
    vector<PseudoJet>::iterator itMother;
    PseudoJet parent1, parent2;
    
    if (fatjet.size() > 0)
    {
        if (fatjet[0].m() > 110)
        {
            mother = {fatjet[0]};
        }
        else
        {
            return {};
        }
    }
    else
    {
        return {};
    }
    
    while (mother.size() > 0)
    {
        for (itMother = mother.begin(); itMother != mother.end(); ++itMother)
        {
            if (!(*itMother).has_parents(parent1, parent2)) continue;
            if (parent1.m() > parent2.m() && parent1.m() > 0.8 * (*itMother).m())
            {
                tmp.push_back(parent1);
            }
            else if (parent2.m() > parent1.m() && parent2.m() > 0.8 * (*itMother).m())
            {
                tmp.push_back(parent2);
            }
            else
            {
                tmp.push_back(parent1);
                tmp.push_back(parent2);
            }
        }
        mother.clear();
        for (int i = 0; i < tmp.size(); i++)
        {
            if (tmp[i].m() < 30)
            {
                subStructure.push_back(tmp[i]);
            }
            else
            {
                mother.push_back(tmp[i]);
            }
            
        }
        tmp.clear();
    }
    //cout << "subs" << " " << subStructure.size() << endl;
    int nJet = subStructure.size();
    double delta;
    vector<PseudoJet> combiedJet, constituent1, constituent2;
    vector<vector<PseudoJet>> higgsCandidate;
    ClusterSequence *filter;
    for (int i = 0; i < nJet; i++)
    {
        for (int j = i + 1; j < nJet; j++)
        {
            constituent1 = subStructure[i].constituents();
            constituent2 = subStructure[j].constituents();
            combiedJet.reserve(constituent1.size() + constituent2.size());
            combiedJet.insert(combiedJet.end(), constituent1.begin(), constituent1.end());
            combiedJet.insert(combiedJet.end(), constituent2.begin(), constituent2.end());
            delta = min(0.3, deltaR(subStructure[i], subStructure[j]));
            filter = new ClusterSequence(combiedJet, JetDefinition(cambridge_algorithm, delta));
            tmp = sorted_by_pt(filter->inclusive_jets());
            if (tmp.size() >= 3)
            {
                higgsCandidate.push_back({tmp[0], tmp[1], tmp[2]});
            }
        }
    }
    return findHiggs(higgsCandidate, parton);
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

    vector<PseudoJet> finalState, higgs;
    vector<GenParticle*> parton;

    int nEvent = treeReader->GetEntries();
    int n = 0;
    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        if(trigger(branchJet))
        {
            finalState = getFinalState(branchTower);
            parton = getParton(branchParticle);
            higgs = higgsTagger(finalState, parton);
            if (higgs.size() == 2)
            {
                if ((higgs[0] + higgs[1]).pt() > 150 && abs((higgs[0] + higgs[1]).m() - 125) < 20)
                {
                    n += 1;
                }
                
            }
        }
        
    }
    cout << n << endl;
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