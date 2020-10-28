#include "BTChannel.h"

using namespace std;
using namespace fastjet;

BTChannel::BTChannel()
{

}

BTChannel::~BTChannel()
{

}

void BTChannel::init(vector<Jet*> jet, vector<Electron*> electron, vector<Muon*> muon, TLorentzVector met)
{
    this->jet = jet;
    this->electron = {};
    this->muon = {};
    this->tauJet = {};
    this->bJet = {};
    this->status = true;
    this->met = met;
    this->lepTau = {};

    for (int i = 0; i < electron.size(); i++)
    {
        if (electron[i]->PT > 20)
        {
            (this->electron).push_back(electron[i]);
        }
    }
    
    for (int i = 0; i < muon.size(); i++)
    {
        if (muon[i]->PT > 18)
        {
            (this->muon).push_back(muon[i]);
        }
    }

    for (int i = 0; i < jet.size(); i++)
    {
        if (jet[i]->TauTag && jet[i]->PT > 20)
        {
            (this->tauJet).push_back(jet[i]);
            jet.erase(jet.begin() + i);
            i--;
            //cout << jet[i]->Mass << endl;
        }
        else if (jet[i]->BTag && jet[i]->PT > 30)
        {
            (this->bJet).push_back(jet[i]);
            //cout << jet[i]->Charge << endl;
            jet.erase(jet.begin() + i);
            i--;
        }
    }

    for (int i = 0; i < jet.size(); i++)
    {
        if (jet[i]->PT > event.hardJet.Pt())
        {
            event.hardJet.SetPtEtaPhiM(jet[i]->PT, jet[i]->Eta, jet[i]->Phi, jet[i]->Mass);
        }
        
    }
    
    //cout << (this->tauJet).size() << endl;
}

bool BTChannel::trigger()
{
    if (tauJet.size() == 0)
    {
        return false;
    }
    else if (tauJet.size() == 1 && bJet.size() == 2)
    {
        if (electron.size() == 1 && muon.size() == 0)
        {
            if (electron[0]->Charge + tauJet[0]->Charge == 0)
            {
                TLorentzVector tmp;
                tmp.SetPtEtaPhiM(electron[0]->PT, electron[0]->Eta, electron[0]->Phi, 0);
                lepTau = {tmp};
                return true;
            }
            else
            {
                return false;
            }
        }
        else if (muon.size() == 1 && electron.size() == 0)
        {
            if (muon[0]->Charge + tauJet[0]->Charge == 0)
            {
                TLorentzVector tmp;
                tmp.SetPtEtaPhiM(muon[0]->PT, muon[0]->Eta, muon[0]->Phi, 0);
                lepTau = {tmp};
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
        
    }
    else if (tauJet.size() == 2 && bJet.size() == 2)
    {
        if ((tauJet[0]->Charge) + (tauJet[1]->Charge) == 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

void BTChannel::findInvMom()
{
    double pNu1, pNu2;
    TLorentzVector pMis1, pMis2, tau1, tau2;
    tau1 = event.tau;
    tau2 = event.antitau;

    if (sin(tau1.Theta()) * sin(tau2.Phi() - tau1.Phi()) == 0 || sin(tau2.Theta()) * sin(tau1.Phi() - tau2.Phi()) == 0)
    {
        event.pMis.SetPxPyPzE(0, 0, 0, 0);
    }
    else
    {
        pNu1 = (met.Px() * sin(tau2.Phi()) - met.Py() * cos(tau2.Phi())) / (sin(tau1.Theta()) * sin(tau2.Phi() - tau1.Phi()));
        pNu2 = (met.Px() * sin(tau1.Phi()) - met.Py() * cos(tau1.Phi())) / (sin(tau2.Theta()) * sin(tau1.Phi() - tau2.Phi()));
        
        pMis1.SetPtEtaPhiM(pNu1 * sin(tau1.Theta()), tau1.Eta(), tau1.Phi(), 0);
        pMis2.SetPtEtaPhiM(pNu2 * sin(tau2.Theta()), tau2.Eta(), tau2.Phi(), 0);

        event.pMis = (pMis1 + pMis2);
    }
}

double BTChannel::mTMax(const double *x)
{
    const double x1 = x[0];
    const double x2 = x[1];
    const double mb = 4.7;
    const double y1 = (met + event.tau + event.antitau).Px() - x1;
    const double y2 = (met + event.tau + event.antitau).Py() - x2;

    const double tmp1 = sqrt(mb * mb + event.tau.M2() + 2 * (sqrt(mb * mb + pow(event.antib.Pt(), 2)) * sqrt(event.tau.M2() + x1 * x1 + x2 * x2) 
                - event.antib.Px() * x1 - event.antib.Py() * x2)); 
    const double tmp2 = sqrt(mb * mb + event.antitau.M2() + 2 * (sqrt(mb * mb + pow(event.b.Pt(), 2)) * sqrt(event.antitau.M2() + y1 * y1 + y2 * y2) 
                - event.b.Px() * y1 - event.b.Py() * y2)); 
    if (tmp1 > tmp2)
    {
        return tmp1;
    }
    else
    {
        return tmp2;
    }
}

void BTChannel::consEvent()
{
    if (tauJet.size() == 2)
    {
        if (tauJet[0]->Charge == 1)
        {
            event.antitau.SetPtEtaPhiM(tauJet[0]->PT, tauJet[0]->Eta, tauJet[0]->Phi, tauJet[0]->Mass);
            event.tau.SetPtEtaPhiM(tauJet[1]->PT, tauJet[1]->Eta, tauJet[1]->Phi, tauJet[1]->Mass);
        }
        else
        {
            event.antitau.SetPtEtaPhiM(tauJet[1]->PT, tauJet[1]->Eta, tauJet[1]->Phi, tauJet[1]->Mass);
            event.tau.SetPtEtaPhiM(tauJet[0]->PT, tauJet[0]->Eta, tauJet[0]->Phi, tauJet[0]->Mass);
        }
    }
    else
    {
        if (tauJet[0]->Charge == 1)
        {
            event.antitau.SetPtEtaPhiM(tauJet[0]->PT, tauJet[0]->Eta, tauJet[0]->Phi, tauJet[0]->Mass);
            event.tau = lepTau[0];
        }
        else
        {
            event.antitau = lepTau[0];
            event.tau.SetPtEtaPhiM(tauJet[0]->PT, tauJet[0]->Eta, tauJet[0]->Phi, tauJet[0]->Mass);
        }
    }
    
    if (bJet[0]->Charge == 1)
    {
        event.antib.SetPtEtaPhiM(bJet[0]->PT, bJet[0]->Eta, bJet[0]->Phi, bJet[0]->Mass);
        event.b.SetPtEtaPhiM(bJet[1]->PT, bJet[1]->Eta, bJet[1]->Phi, bJet[1]->Mass);
    }
    else
    {
        event.b.SetPtEtaPhiM(bJet[0]->PT, bJet[0]->Eta, bJet[0]->Phi, bJet[0]->Mass);
        event.antib.SetPtEtaPhiM(bJet[1]->PT, bJet[1]->Eta, bJet[1]->Phi, bJet[1]->Mass);
    }
    findInvMom();
}

bool BTChannel::selector()
{
    double mtt = (event.tau + event.antitau + event.pMis).M();
    double mbb = (event.b + event.antib).M();
    double deltaBB = event.b.DeltaR(event.antib);
    double ptBB = (event.b + event.antib).Pt();
    double ptj = event.hardJet.Pt();

    if (mtt > 80 && mtt < 170 && mbb > 100 && mbb < 150 && deltaBB > 0.4 && deltaBB < 2 /* && ptBB > 80 */ && ptj > 120)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void BTChannel::process()
{
    if (trigger())
    {
        consEvent();
        if (selector())
        {
            dihiggsInvM = (event.tau + event.antitau + event.b + event.antib + event.pMis).M();
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