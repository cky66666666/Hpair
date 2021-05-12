#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TF1.h"

#include <stdlib.h>
#include "iostream"
#include "vector"
#include "fstream"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

using namespace std;

vector<TLorentzVector> ptSort(vector<TLorentzVector> lightJet)
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

vector<TLorentzVector> wReco(vector<TLorentzVector> lightJet)
{
    if (lightJet.size() < 2) return {};
    double deltaM = abs((lightJet[0] + lightJet[1]).M() - 80);
    vector<TLorentzVector> wCand = {lightJet[0], lightJet[1]};
    for (int i = 0; i < lightJet.size(); i++)
    {
        for (int j = i + 1; j < lightJet.size(); j++)
        {
            if (abs((lightJet[i] + lightJet[j]).M() - 80) < deltaM)
            {
                deltaM = abs((lightJet[i] + lightJet[j]).M() - 80);
                wCand = {lightJet[i], lightJet[j]};
            }
            
        }
        
    }
    return wCand;
}

vector<TLorentzVector> topReco(vector<TLorentzVector> bJet, vector<TLorentzVector> lightJet)
{
    if (bJet.size() < 1 || lightJet.size() < 2) return {};
    vector<TLorentzVector> wCand = wReco(lightJet);
    TLorentzVector b = bJet[0];
    double deltaM = abs((b + wCand[0] + wCand[1]).M() - 173);

    for (int i = 0; i < bJet.size(); i++)
    {
        if (abs((bJet[i] + wCand[0] + wCand[1]).M() - 173) <  deltaM)
        {
            deltaM = abs((bJet[i] + wCand[0] + wCand[1]).M() - 173);
            b = bJet[i];
        }
    }
    wCand.push_back(b);
    return wCand;
}

int main()
{
    TChain *chain = new TChain("Delphes");
    chain->Add("/home/E/chaikangyu/work/Hpair/events/root/bkg_aa_tthj_100_100_1.root");
    chain->Add("/home/E/chaikangyu/work/Hpair/events/root/bkg_aa_tthj_100_100_2.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");

    TH1D *histw = new TH1D("w", "w", 50, 30, 130);
    TH1D *histtop = new TH1D("t", "t", 50, 120, 220);
    Jet *jet;
    for (int i = 0; i < treeReader->GetEntries(); i++)
    {
        treeReader->ReadEntry(i);
        vector<TLorentzVector> bJet = {}, lightJet = {};
        TLorentzVector tmp;
        vector<TLorentzVector> topCand;
        for (int j = 0; j < branchJet->GetEntries(); j++)
        {
            jet = (Jet*) branchJet->At(j);
            tmp.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
            if (jet->BTag == 1 && jet->PT > 30)
            {
                bJet.push_back(tmp);
            }
            else if (jet->BTag != 1 && jet->PT > 30)
            {
                lightJet.push_back(tmp);
            }
        }
        /* topCand = topReco(bJet, lightJet);
        if(topCand.size() == 3)
        {
            histw->Fill((topCand[0] + topCand[1]).M());
            histtop->Fill((topCand[0] + topCand[1] + topCand[2]).M());
        } */
        lightJet = ptSort(lightJet);
        bJet = ptSort(bJet);
        if (lightJet.size() >=3 && bJet.size() >= 2)
        {
            histtop->Fill((lightJet[1] + lightJet[2] +bJet[0]).M());
            histw->Fill((lightJet[1] + lightJet[2]).M());
        }
        
    }
    TF1 *gaus = new TF1("gaus", "gaus", 30, 130);
    histw->Fit(gaus);
    histtop->Fit(gaus);
    TFile *f1 = new TFile("w.root", "RECREATE");
    histw->Write();
    f1->Close();
    TFile *f2 = new TFile("t.root", "RECREATE");
    histtop->Write();
    f2->Close();
    return 0;
}