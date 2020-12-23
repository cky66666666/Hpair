#include "BAChannel.h"

using namespace std;
using namespace fastjet;


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

void BAChannel::init(vector<PseudoJet> finalState, vector<GenParticle*> parton, vector<TLorentzVector> photon, vector<TLorentzVector> electron, vector<TLorentzVector> muon, vector<Jet*> delphesJet)
{
    this->finalState = finalState;
    this->parton = parton;
    this->photon = photon;
    this->electron = electron;
    this->muon = muon;
    status = true;
    higgsCandListA = {};
    higgsCandListB = {};
    bJet = {};
    lightJet = {};
    this->delphesJet = delphesJet;
}

bool BAChannel::trigger()
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
    TLorentzVector pJet;
    
    /* ClusterSequence sequence = ClusterSequence(finalState, JetDefinition(antikt_algorithm, 0.4));
    vector<PseudoJet> tmpJet = sorted_by_pt(sequence.inclusive_jets(30));

    for (int i = 0; i < tmpJet.size(); i++)
    {
        if (tmpJet[i].pt() > 30 && abs(tmpJet[i].eta()) < 2.5)
        {
            pJet.SetPxPyPzE(tmpJet[i].px(), tmpJet[i].py(), tmpJet[i].pz(), tmpJet[i].e());
            if (flavourAssociation(pJet) == 5)
            {
                bJet.push_back(pJet);
            }
            else
            {
                lightJet.push_back(pJet);
            }
            
        }
    } */
    double Ht = 0;
    for (int i = 0; i < delphesJet.size(); i++)
    {
        if (delphesJet[i]->PT > 30 && abs(delphesJet[i]->Eta) < 2.5)
        {
            pJet.SetPtEtaPhiM(delphesJet[i]->PT, delphesJet[i]->Eta, delphesJet[i]->Phi, delphesJet[i]->Mass);
            Ht += pJet.Pt();
            if (delphesJet[i]->BTag == 1)
            {
                bJet.push_back(pJet);
            }
            else
            {
                lightJet.push_back(pJet);
            }
        }
        
    }
    signal.nJet = bJet.size() + lightJet.size();
    signal.Ht = Ht;
    
    for (int i = 0; i < photon.size(); i++)
    {
        if (photon[i].Pt() < 30 || (abs(photon[i].Eta()) < 1.52 && abs(photon[i].Eta()) > 1.37) || abs(photon[i].Eta()) >= 2.5)
        {
            photon.erase(photon.begin() + i);
            i--;
        }
    }

}

void BAChannel::selPhotonPair()
{
    int nPhoton = photon.size();
    vector<double> diPhotonPt = {};
    vector<double> deltaInvM = {};

    for (int i = 0; i < nPhoton; i++)
    {
        for (int j = i + 1; j < nPhoton; j++)
        {
            if ((photon[i].DeltaR(photon[j]) > 0.4) && (photon[i].DeltaR(photon[j]) < 2.0) && ((photon[i] + photon[j]).M() < 128) && ((photon[i] + photon[j]).M() > 122))
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
    
    if (higgsCandListA.size() > 0)
    {
        for (int i = 0; i < higgsCandListA.size(); i++)
        {
            diPhotonPt.push_back((higgsCandListA[i].cand1 + higgsCandListA[i].cand2).Pt());
            deltaInvM.push_back(abs((higgsCandListA[i].cand1 + higgsCandListA[i].cand2).M() - 125));
        }
        int maxPosition = max_element(diPhotonPt.begin(), diPhotonPt.end()) - diPhotonPt.begin();
        photonPair = {higgsCandListA[maxPosition].cand1, higgsCandListA[maxPosition].cand2};
    }
    

}

void BAChannel::selBPair()
{
    vector<double> bPairPt = {};
    vector<double> deltaInvM = {};
    for (int i = 0; i < bJet.size(); i++)
    {
        for (int j = i + 1; j < bJet.size(); j++)
        {
            if ((bJet[i].DeltaR(bJet[j]) < 2) && (bJet[i].DeltaR(bJet[j]) > 0.4) 
                && ((bJet[i] + bJet[j]).M() < 150) && ((bJet[i] + bJet[j]).M() > 100) && max(bJet[i].Pt(), bJet[j].Pt()) > 40 && bJet[i].DeltaR(photonPair[0]) > 0.4
                && bJet[i].DeltaR(photonPair[1]) > 0.4 && bJet[j].DeltaR(photonPair[0]) > 0.4 && bJet[j].DeltaR(photonPair[1]) > 0.4)
            {
                HiggsCand tmp;
                tmp.cand1 = bJet[i];
                tmp.cand2 = bJet[j];
                tmp.index1 = i;
                tmp.index2 = j;
                higgsCandListB.push_back(tmp);
            }
        }   
    }

    if (higgsCandListB.size() > 0)
    {
        for (int i = 0; i < higgsCandListB.size(); i++)
        {
            bPairPt.push_back((higgsCandListB[i].cand1 + higgsCandListB[i].cand2).Pt());
            deltaInvM.push_back(abs((higgsCandListB[i].cand1 + higgsCandListB[i].cand2).M() - 125));
        }
        int maxPosition = max_element(bPairPt.begin(), bPairPt.end()) - bPairPt.begin();
        bPair = {higgsCandListB[maxPosition].cand1, higgsCandListB[maxPosition].cand2};
        bJet.erase(bJet.begin() + higgsCandListB[maxPosition].index2);
        bJet.erase(bJet.begin() + higgsCandListB[maxPosition].index1);
    }

    //vector<TLorentzVector> tmp;
    vector<vector<TLorentzVector>> tmp1 = {bJet, lightJet};
    lightJet = combineVector(tmp1);
}

void BAChannel::find2BHiggsHard()
{
    BoostedHiggs boostedHiggs;
    vector<PseudoJet> remmant, remmantJet, photonPair_J;

    for (int i = 0; i < photonPair.size(); i++)
    {
        PseudoJet tmp = PseudoJet(photonPair[0].Px(), photonPair[0].Py(), photonPair[0].Pz(), photonPair[0].E());
        photonPair_J.push_back(tmp);
    }

    boostedHiggs.init(finalState, parton, photonPair_J, 100, 100, 1.5, 0.3);
    boostedHiggs.process(1);
    if (boostedHiggs.boostedHiggs.size() != 0)
    {
        signal.higgs1.SetPxPyPzE((boostedHiggs.boostedHiggs)[0].px(), (boostedHiggs.boostedHiggs)[0].py(), (boostedHiggs.boostedHiggs)[0].pz(), (boostedHiggs.boostedHiggs)[0].e());
        remmant = boostedHiggs.remmant;
        ClusterSequence remmantSeq = ClusterSequence(remmant, JetDefinition(antikt_algorithm, 0.4));
        remmantJet = sorted_by_pt(remmantSeq.inclusive_jets(30));
        if (remmantJet.size() != 0)
        {
            signal.hardJet.SetPxPyPzE(remmantJet[0].px(), remmantJet[0].py(), remmantJet[0].pz(), remmantJet[0].e());
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
    //HiggsCand higgsA, higgsB;
    preprocess();
    if (trigger())
    {
        selPhotonPair();
        if (photonPair.size() == 2)
        {
            selBPair();
            if (bPair.size() == 2 && lightJet.size() > 0)
            {
                signal.b1 = bPair[0];
                signal.b2 = bPair[1];
                signal.other1 = photonPair[0];
                signal.other2 = photonPair[1];
                signal.higgs1 = bPair[0] + bPair[1];
                signal.higgs2 = photonPair[0] + photonPair[1];

                vector<double> jetPt = {};
                for (int i = 0; i < lightJet.size(); i++)
                {
                    jetPt.push_back(lightJet[i].Pt());
                }
                int maxPosition = max_element(jetPt.begin(), jetPt.end()) - jetPt.begin();
                //cout << jetPt[maxPosition] << endl;
                signal.hardJet = lightJet[maxPosition];
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
    bJet.clear();
    lightJet.clear();
    higgsCandListA.clear();
    higgsCandListB.clear();
}