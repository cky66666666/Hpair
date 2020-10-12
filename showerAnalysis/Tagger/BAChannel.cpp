#include "BAChannel.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "classes/DelphesClasses.h"
#include "vector"

using namespace std;
using namespace fastjet;
struct HiggsCand
{
    TLorentzVector cand1;
    TLorentzVector cand2;
    int index1, index2;
};

int BAChannel::flavourAssociation(TLorentzVector jet)
{
    int flavour = -1;
    double deltaR = 0.5;
    for (int i = 0; i < parton.size(); i++)
    {
        int pid = abs(parton[i]->PID);
        if (pid == 21) pid = 0;
        if (jet.DeltaR(parton[i]->P4()) < deltaR && pid > flavour)
        {
            flavour = pid;
        }
    }
    return flavour;
}

void BAChannel::init(vector<PseudoJet> finalState, vector<GenParticle*> parton, vector<TLorentzVector> photon, vector<TLorentzVector> electron, vector<TLorentzVector> muon)
{
    this->finalState = finalState;
    this->parton = parton;
    this->photon = photon;
    this->electron = electron;
    this->muon = muon;
    status = true;
    higgsCandListA = {};
    higgsCandListB = {};
    jet = {};
    bJet = {};
}

void BAChannel::preprocess()
{
    ClusterSequence sequence = ClusterSequence(finalState, JetDefinition(antikt_algorithm, 0.4));
    vector<PseudoJet> tmpJet = sorted_by_pt(sequence.inclusive_jets(30));
    TLorentzVector pJet;

    for (int i = 0; i < jet.size(); i++)
    {
        if (tmpJet[i].pt() > 30 && abs(tmpJet[i].eta()) < 2.5)
        {
            pJet.SetPxPyPzE(tmpJet[i].px(), tmpJet[i].py(), tmpJet[i].pz(), tmpJet[i].e());
            jet.push_back(pJet);
        }
    }

    for (int i = 0; i < photon.size(); i++)
    {
        if (photon[i].Pt() > 30 || !(abs(photon[i].Pt()) < 1.37 || 1.52 < abs(photon[i].Pt()) < 2.37))
        {
            photon.erase(photon.begin() + i);
            i--;
        }
    }

}

void BAChannel::find2AHiggs()
{
    int nPhoton = photon.size();

    for (int i = 0; i < nPhoton; i++)
    {
        for (int j = i + 1; j < nPhoton; j++)
        {
            if ((0.4 < photon[i].DeltaR(photon[j]) < 2.0) && (122 < (photon[i] + photon[j]).M() < 128))
            {
                HiggsCand tmp;
                tmp.cand1 = photon[i];
                tmp.cand2 = photon[j];
                higgsCandListA.push_back(tmp);
            }
        } 
    }
}

void BAChannel::find2BHiggs()
{
    for (int i = 0; i < jet.size(); i++)
    {
        for (int j = i + 1; j < jet.size(); j++)
        {
            if (flavourAssociation(jet[i]) == 5 && flavourAssociation(jet[j]) == 5 && (0.4 < jet[i].DeltaR(jet[j]) < 2) && (100 < (jet[i] + jet[j]).M() < 150) && jet[i].Pt() > 40)
            {
                HiggsCand tmp;
                tmp.cand1 = jet[i];
                tmp.cand2 = jet[j];
                tmp.index1 = i;
                tmp.index2 = j;
                higgsCandListB.push_back(tmp);
            }
            
        }   
    }
}