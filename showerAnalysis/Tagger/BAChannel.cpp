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

void BAChannel::init(vector<PseudoJet> finalState, vector<GenParticle*> parton, vector<TLorentzVector> photon, 
                vector<TLorentzVector> electron, vector<TLorentzVector> muon, vector<Jet*> delphesJet, Cut cut)
{
    this->finalState = finalState;
    this->parton = parton;
    this->photon = photon;
    this->electron = electron;
    this->muon = muon;
    this->cut = cut;
    status = true;
    higgsCandListA = {};
    higgsCandListB = {};
    bJet = {};
    lightJet = {};
    bPair = {};
    photonPair ={};
    this->delphesJet = delphesJet;
}

double BAChannel::thrustTarget(const double *angle)
{
    const double phi = angle[0];
    const double theta = angle[1];
    TVector3 nHat(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    double result = 0, norm = 0;
    vector<vector<TLorentzVector>> tmp = {photon, electron, muon, lightJet, bJet};
    vector<TLorentzVector> objectList = combineVector(tmp);
    for (int i = 0; i < objectList.size(); i++)
    {
        result += objectList[i].Vect().Dot(nHat);
        norm += objectList[i].Vect().Mag();
    }
    return -result / norm;
}

double BAChannel::chi2(TLorentzVector j1, TLorentzVector j2, TLorentzVector b)
{
    return pow((j1 + j2).M() - 81.8653, 2) / 10.81 / 10.81 + pow((j1 + j2 + b).M() - 181.505, 2) / 31.0125 / 31.0125;
}

double BAChannel::topness()
{
    double top = chi2(lightJet[1], lightJet[2], bPair[0]);
    
    for (int i = 0; i < bPair.size(); i++)
    {
        for (int j = 1; j < lightJet.size(); j++)
        {
            for (int k = j + 1; k < lightJet.size(); k++)
            {
                double tmp = chi2(lightJet[j], lightJet[k], bPair[i]);
                if (tmp < top)
                {
                    top = tmp;
                }
            }
        }
    }
    return top;
}

double BAChannel::topness2()
{
    double top = chi2(lightJet[1], lightJet[2], bPair[0]) + chi2(lightJet[3], lightJet[4], bPair[1]);

    for (int i = 0; i < bPair.size(); i++)
    {
        for (int j = 1; j < lightJet.size(); j++)
        {
            for (int k = j + 1; k < lightJet.size(); k++)
            {
                vector<TLorentzVector> tmp1 = lightJet;
                tmp1.erase(tmp1.begin() + k);
                tmp1.erase(tmp1.begin() + j);
                for (int l = 1; l < tmp1.size(); l++)
                {
                    for (int m = l + 1; m < tmp1.size(); m++)
                    {
                        double tmp;
                        if (i == 0)
                        {
                            tmp = chi2(lightJet[j], lightJet[k], bPair[0]) + chi2(tmp1[l], tmp1[m], bPair[1]);
                        }
                        else
                        {
                            tmp = chi2(lightJet[j], lightJet[k], bPair[1]) + chi2(tmp1[l], tmp1[m], bPair[0]);
                        }
                        if (tmp < top)
                        {
                            top = tmp;
                        }
                    }
                }
            }
        }
    }
    return top;
}

double BAChannel::dihedralAngle(vector<TLorentzVector> surf1, vector<TLorentzVector> surf2)
{
    if (surf1.size() != 2 || surf2.size() != 2) return 0;
    TVector3 v1 = surf1[0].Vect(), v2 = surf1[1].Vect(), w1 = surf2[0].Vect(), w2 = surf2[1].Vect();
    TVector3 n1 = v1.Cross(v2), n2 = w1.Cross(w2);
    return n1.Angle(n2);
}

double BAChannel::deltaPhi(vector<TLorentzVector> surf, TLorentzVector vec)
{
    if (surf.size() != 2) return 0;
    TVector3 v1 = surf[0].Vect(), v2 = surf[1].Vect(), w = vec.Vect();
    return (v1.Cross(v2)).Angle(w);
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
            if ((photon[i].DeltaR(photon[j]) > cut.deltaR_aa_min) && (photon[i].DeltaR(photon[j]) < cut.deltaR_aa_max) 
                && abs((photon[i] + photon[j]).M() - 125) < cut.delta_maa)
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
            if ((bJet[i].DeltaR(bJet[j]) < cut.deltaR_bb_max) && (bJet[i].DeltaR(bJet[j]) > cut.deltaR_bb_min) 
                && abs((bJet[i] + bJet[j]).M() - 125) < cut.delta_mbb && max(bJet[i].Pt(), bJet[j].Pt()) > cut.ptb 
                && bJet[i].DeltaR(photonPair[0]) > cut.deltaR_ab_min && bJet[i].DeltaR(photonPair[1]) > cut.deltaR_ab_min 
                && bJet[j].DeltaR(photonPair[0]) > cut.deltaR_ab_min && bJet[j].DeltaR(photonPair[1]) > cut.deltaR_ab_min)
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

vector<TLorentzVector> BAChannel::ptSort(vector<TLorentzVector> lightJet)
{
    vector<TLorentzVector> result = lightJet;
    int n = result.size();

    while (n > 1)
    {
        for (int i = 0; i < n - 1; i++)
        {
            if (result[i].Pt() < result[i + 1].Pt())
            {
                TLorentzVector tmp = result[i];
                result[i] = result[i + 1];
                result[i + 1] = tmp;
            }
            n--;
        }
    }
    return result;
}

void BAChannel::process()
{
    //HiggsCand higgsA, higgsB;
    preprocess();
    signal.status = true;
    if (trigger())
    {
        selPhotonPair();
        if (photonPair.size() == 2)
        {
            selBPair();
        }        
        if (photonPair.size() == 2 && bPair.size() == 2 && lightJet.size() > 0)
        {
            if (/* (bPair[0].Pt() > bPair[1].Pt() && bPair[0].Pt() > 40 && bPair[1].Pt() > 30) 
            || (bPair[0].Pt() < bPair[1].Pt() && bPair[0].Pt() > 30 && bPair[1].Pt() > 40) */ true)
            {
                bPair = ptSort(bPair);
                lightJet = ptSort(lightJet);
                signal.b1 = bPair[0];
                signal.b2 = bPair[1];
                signal.other1 = photonPair[0];
                signal.other2 = photonPair[1];
                signal.higgs1 = bPair[0] + bPair[1];
                signal.higgs2 = photonPair[0] + photonPair[1];

                signal.jetList = lightJet;
                signal.hardJet = lightJet[0];
                signal.dih_bbaa = dihedralAngle(bPair, photonPair);
                signal.deltaPhi = (photonPair[0] + photonPair[1]).DeltaPhi(bPair[1]);
                signal.nJet = bPair.size() + lightJet.size();
                if (lightJet.size() < 3)
                {
                    signal.topness = -1;
                    signal.topness2 = -1;
                }
                else if (lightJet.size() >= 3 && lightJet.size() < 5)
                {
                    signal.topness = topness();
                    signal.topness2 = -1;
                }
                else
                {
                    signal.topness = topness();
                    signal.topness2 = topness2();
                }
            }
            else
            {
                signal.status = false;
            }
        }
        else
        {
            signal.status = false;
        }
    }
    else
    {
        signal.status = false;
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