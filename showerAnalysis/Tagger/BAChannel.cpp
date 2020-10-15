#include "BAChannel.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "classes/DelphesClasses.h"
#include "vector"
#include "BoostedHiggs.h"

using namespace std;
using namespace fastjet;

class HiggsCand
{
public:
    HiggsCand()
    {

    }
    ~HiggsCand()
    {

    }
    TLorentzVector cand1;
    TLorentzVector cand2;
    int index1, index2;
};

BAChannel::BAChannel()
{

}

BAChannel::~BAChannel()
{

}

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

HiggsCand BAChannel::higgsSelector(vector<HiggsCand> candidate)
{
    int n = candidate.size();
    HiggsCand higgsCand;
    
    if (n == 1)
    {
        higgsCand = candidate[0];
    }
    else
    {
        TLorentzVector tmp;
        double massDif = 1000;
        for (int i = 0; i < n; i++)
        {
            tmp = candidate[i].cand1 + candidate[i].cand2;
            if (abs(tmp.M() - 125) < massDif)
            {
                higgsCand = candidate[i];
            }
        }
    }
    if (higgsCand.index1 != higgsCand.index2)
    {
        if (higgsCand.index1 > higgsCand.index2)
        {
            jet.erase(jet.begin() + higgsCand.index1);
            jet.erase(jet.begin() + higgsCand.index2);
        }
        else
        {
            jet.erase(jet.begin() + higgsCand.index2);
            jet.erase(jet.begin() + higgsCand.index1);
        }
    }
    return higgsCand;
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
}

bool BAChannel::preselect()
{
    for (int i = 0; i < electron.size(); i++)
    {
        if (electron[i].Pt() > 25 && abs(electron[i].Eta()) < 2.5)
        {
            return false;
        }
    }

    for (int i = 0; i < muon.size(); i++)
    {
        if (muon[i].Pt() > 25 && abs(muon[i].Eta()) < 2.5)
        {
            return false;
        }
    }
    
    if (photon.size() < 2 /* || jet.size() < 3 || jet.size() > 7 */)
    {
        return false;
    }
    else
    {
        return true;
    }
    
    
}

void BAChannel::preprocess()
{
    /* ClusterSequence sequence = ClusterSequence(finalState, JetDefinition(antikt_algorithm, 0.4));
    vector<PseudoJet> tmpJet = sorted_by_pt(sequence.inclusive_jets(30));
    TLorentzVector pJet; */

    /* for (int i = 0; i < tmpJet.size(); i++)
    {
        if (tmpJet[i].pt() > 30 && abs(tmpJet[i].eta()) < 2.5)
        {
            pJet.SetPxPyPzE(tmpJet[i].px(), tmpJet[i].py(), tmpJet[i].pz(), tmpJet[i].e());
            jet.push_back(pJet);
        }
    } */

    for (int i = 0; i < photon.size(); i++)
    {
        if (photon[i].Pt() < 30 /* || !(abs(photon[i].Eta()) < 1.37 || 1.52 < abs(photon[i].Eta()) < 2.37) */)
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
            if ((0.4 < photon[i].DeltaR(photon[j]) < 2.0) && (123 < (photon[i] + photon[j]).M() < 127))
            {
                HiggsCand tmp;
                tmp.cand1 = photon[i];
                tmp.cand2 = photon[j];
                tmp.index1 = 0;
                tmp.index2 = 0;
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

void BAChannel::find2BHiggsHard()
{
    BoostedHiggs boostedHiggs;
    vector<PseudoJet> remmant, remmantJet;
    
    boostedHiggs.init(finalState, parton, 100, 110, 2, 0.3);
    boostedHiggs.process(1);
    if (boostedHiggs.boostedHiggs.size() != 0)
    {
        higgsFromB.SetPxPyPzE((boostedHiggs.boostedHiggs)[0].px(), (boostedHiggs.boostedHiggs)[0].py(), (boostedHiggs.boostedHiggs)[0].pz(), (boostedHiggs.boostedHiggs)[0].e());
        remmant = boostedHiggs.remmant;
        ClusterSequence remmantSeq = ClusterSequence(remmant, JetDefinition(antikt_algorithm, 0.4));
        remmantJet = sorted_by_pt(remmantSeq.inclusive_jets(30));
        if (remmantJet.size() != 0)
        {
            hardJet.SetPxPyPzE(remmantJet[0].px(), remmantJet[0].py(), remmantJet[0].pz(), remmantJet[0].e());
        }
        else
        {
            status = false;
        }
    }
    else
    {
        status = false;
    }
    boostedHiggs.clear();

}

void BAChannel::process()
{
    HiggsCand higgsA, higgsB;
    preprocess();
    if (preselect())
    {
        find2AHiggs();
        find2BHiggsHard();

        /* if (higgsCandListA.size() == 0 || higgsCandListB.size() == 0)
        {
            status = false;
        }
        else
        {
            higgsA = higgsSelector(higgsCandListA);
            higgsB = higgsSelector(higgsCandListB);

            if (higgsA.cand1.DeltaR(higgsB.cand1) < 0.4 || higgsA.cand1.DeltaR(higgsB.cand2) < 0.4 
                || higgsA.cand2.DeltaR(higgsB.cand1) < 0.4 || higgsA.cand2.DeltaR(higgsB.cand2) < 0.4)
            {
                status = false;
            }
            else
            {
                hardJet = jet[0];
                higgsFromA = higgsA.cand1 + higgsA.cand2;
                higgsFromB = higgsB.cand1 + higgsB.cand2;
            }
        } */
        if (higgsCandListA.size() == 0)
        {
            status = false;
        }
        else
        {
            higgsA = higgsSelector(higgsCandListA);
            higgsFromA = higgsA.cand1 + higgsA.cand2;
        }
        
        
    }
    else
    {
        status = false;
    }
}

void BAChannel::finish()
{
    finalState.clear();
    parton.clear();
    photon.clear();
    electron.clear();
    muon.clear();
    jet.clear();
    higgsCandListA.clear();
    higgsCandListB.clear();
}