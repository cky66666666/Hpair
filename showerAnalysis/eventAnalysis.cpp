#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TObject.h"

#include <stdlib.h>
#include "iostream"
#include "vector"
#include "fstream"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "Tagger/BoostedHiggs.h"
#include "Tagger/BAChannel.h"
#include "Tagger/BTChannel.h"
#include "Tagger/BasicObject.h"

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
        jetPtVec.push_back(jet->PT);
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

vector<TLorentzVector> getObject(TClonesArray *branchObject, int objType)
{
    //objType=1:photon
    //objType=2:electron
    //objType=3:muon

    int nObj = branchObject->GetEntries();
    vector<TLorentzVector> object = {};
    TLorentzVector momentum;

    if (objType == 1)
    {
        Photon *photon;
        for (int i = 0; i < nObj; i++)
        {
            photon = (Photon*) branchObject->At(i);
            momentum.SetPtEtaPhiE(photon->PT, photon->Eta, photon->Phi, photon->E);
            object.push_back(momentum);
        }
        return object;
    }
    else if (objType == 2)
    {
        Electron *electron;
        for (int i = 0; i < nObj; i++)
        {
            electron = (Electron*) branchObject->At(i);
            momentum.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0);
            object.push_back(momentum);
        }
        return object;
    }
    else if (objType == 3)
    {
        Muon *muon;
        for (int i = 0; i < nObj; i++)
        {
            muon = (Muon*) branchObject->At(i);
            momentum.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0);
            object.push_back(momentum);
        }
        return object;
    }
    else
    {
        cout << "error type should be 1 2 or 3" << endl;
        return {};
    }
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

void jetFakePhoton(vector<Jet*> &delphesJet, vector<TLorentzVector> &photon, int nFake)
{
    vector<Jet*> lightJet, heavyJet;
    TLorentzVector jetPhoton;

    for (int i = 0; i < delphesJet.size(); i++)
    {
        if (delphesJet[i]->BTag == 0)
        {
            lightJet.push_back(delphesJet[i]);
        }
        else
        {
            heavyJet.push_back(delphesJet[i]);
        }
    }
    if (lightJet.size() < nFake)
    {
        return;
    }
    
    for (int i = 0; i < nFake; i++)
    {
        int n = rand() % (lightJet.size());
        jetPhoton.SetPtEtaPhiM(lightJet[n]->PT, lightJet[n]->Eta, lightJet[n]->Phi, lightJet[n]->Mass);
        photon.push_back(jetPhoton);
        lightJet.erase(lightJet.begin() + n);
    }
    
    vector<Jet*> tmp;
    tmp.reserve(lightJet.size() + heavyJet.size());
    tmp.insert(tmp.end(), lightJet.begin(), lightJet.end());
    tmp.insert(tmp.end(), heavyJet.begin(), heavyJet.end());
    delphesJet = tmp;
}

void fakeBJet(vector<Jet*> delphesJet, int nFake)
{
    int nJet = 0;
    for (int i = 0; i < delphesJet.size(); i++)
    {
        if (delphesJet[i]->BTag != 1)
        {
            nJet += 1;
        }
    }
    if (nJet < nFake) return;
    
    int fake = 0;
    while (fake < nFake)
    {
        int iFake = rand() % (delphesJet.size());
        if (delphesJet[iFake]->BTag != 1)
        {
            delphesJet[iFake]->BTag = 1;
            fake += 1;
        }
    }
    
}

SignalEvent analyseAA(TClonesArray *branchParticle, TClonesArray *branchTower, TClonesArray *branchPhoton, TClonesArray *branchElectron, 
            TClonesArray *branchMuon, TClonesArray *branchJet, Cut cut)
{
    BAChannel *BAEvent = new BAChannel();
    SignalEvent event;
    vector<PseudoJet> finalState;
    vector<GenParticle*> parton;
    vector<TLorentzVector> photon, electron, muon;
    vector<Jet*> delphesJet = {};
    double inv;

    finalState = getFinalState(branchTower);
    parton = getParton(branchParticle);
    photon = getObject(branchPhoton, 1);
    electron = getObject(branchElectron, 2);
    muon = getObject(branchMuon, 3);

    for (int i = 0; i < branchJet->GetEntries(); i++)
    {
        delphesJet.push_back((Jet*) branchJet->At(i));
    }
    
    jetFakePhoton(delphesJet, photon, cut.fakePhoton);
    fakeBJet(delphesJet, cut.fakebJet);

    BAEvent->init(finalState, parton, photon, electron, muon, delphesJet, cut);
    BAEvent->process();
    event = BAEvent->signal;
    /* ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer();
    ROOT::Math::Functor f(BAEvent, &BAChannel::thrustTarget, 2);
    
    minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimizer->SetMaxIterations(10000);  // for GSL
    minimizer->SetTolerance(0.001);
    minimizer->SetPrintLevel(0);
    minimizer->SetFunction(f);

    minimizer->SetVariable(0, "phi", 1, 0.01);
    minimizer->SetVariable(1, "theta", 1, 0.01);
    minimizer->Minimize();

    event.thrust = -minimizer->MinValue(); */
    BAEvent->finish();
    delete BAEvent;
    return event;
}

int main(int argc, char *argv[])
{
    // Usage: ./eventAnalysis /path/to/hist.root /path/to/inputfile numberOfInputFile histogram name cutFile

    //gSystem->Load("/mnt/d/work/Hpair/Delphes/libDelphes");

    /* TH1D *hist = new TH1D(argv[argc - 2], argv[argc - 2], 30, 250, 1000);
    TFile *f = new TFile(argv[1], "RECREATE"); */

    TH1D *histInvM = new TH1D((string(argv[argc - 2]) + string("_InvMass")).c_str(), "InvMass", 50, 250, 1750);
    TH1D *histDeltaPhi = new TH1D((string(argv[argc - 2]) + string("_DeltaPhi")).c_str(), "DeltaPhi", 30, 0, pi);
    TH1D *histBJetRatio = new TH1D((string(argv[argc - 2]) + string("_BJetRatio")).c_str(), "BJetRatio", 20, 0, 1);
    TH1D *histNLJet = new TH1D((string(argv[argc - 2]) + string("_NLJet")).c_str(), "NLJet", 20, 0, 1);
    TH1D *histN2LJet = new TH1D((string(argv[argc - 2]) + string("_N2LJet")).c_str(), "N2LJet", 20, 0, 1);
    TH1D *histN3LJet = new TH1D((string(argv[argc - 2]) + string("_N3LJet")).c_str(), "N3LJet", 20, 0, 1);
    TH1D *histTopness = new TH1D((string(argv[argc - 2]) + string("_Topness")).c_str(), "Topness", 100, 0, 50);
    TH1D *histTopness2 = new TH1D((string(argv[argc - 2]) + string("_Topness2")).c_str(), "Topness2", 100, 0, 50);
    TH1D *histThrust = new TH1D((string(argv[argc - 2]) + string("_Thrust")).c_str(), "Thrust", 20, 0, 1);
    TH1D *histDihedralAngle = new TH1D((string(argv[argc - 2]) + string("_DihedralAngle")).c_str(), "DihedralAngle", 20, 0, 1);
    TH1D *histHRatio = new TH1D((string(argv[argc - 2]) + string("_HRatio")).c_str(), "HRatio", 20, 0, 4);

    TChain *chain = new TChain("Delphes");
    int nData = (int) *argv[3] - 48;

    if (nData == 1)
    {
        chain->Add(argv[2]);
    }
    else
    {
        for (int i = 1; i < nData + 1; i++)
        {
            char inputFile[100];
            sprintf(inputFile, "%s_%d.root", argv[2], i);
            chain->Add(inputFile);
        }
    }

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");

    int nEvent = treeReader->GetEntries();
    //cout << nEvent << endl;
    Cut cut;

    ifstream cutFile(argv[argc - 1]);
    string tmp, cutName, cutValue;
    SignalEvent event;

    while (getline(cutFile, tmp))
    {
        int splitIndex = tmp.find_first_of(" ");
        cutName = tmp.substr(0, splitIndex);
        cutValue = tmp.substr(splitIndex + 1, tmp.length());
        if (cutName == "deltaR_bb_max")
        {
            cut.deltaR_bb_max = atof(cutValue.c_str());
        }
        else if (cutName == "deltaR_bb_min")
        {
            cut.deltaR_bb_min = atof(cutValue.c_str());
        }
        else if (cutName == "deltaR_aa_max")
        {
            cut.deltaR_aa_max = atof(cutValue.c_str());
        }
        else if (cutName == "deltaR_aa_min")
        {
            cut.deltaR_aa_min = atof(cutValue.c_str());
        }
        else if (cutName == "deltaR_ab_max")
        {
            cut.deltaR_ab_max = atof(cutValue.c_str());
        }
        else if (cutName == "deltaR_ab_min")
        {
            cut.deltaR_ab_min = atof(cutValue.c_str());
        }
        else if (cutName == "delta_maa")
        {
            cut.delta_maa = atof(cutValue.c_str());
        }
        else if (cutName == "delta_mbb")
        {
            cut.delta_mbb = atof(cutValue.c_str());
        }
        else if (cutName == "ptb")
        {
            cut.ptb = atof(cutValue.c_str());
        }
        else if (cutName == "ptj")
        {
            cut.ptj = atof(cutValue.c_str());
        }
        else if (cutName == "fakePhoton")
        {
            cut.fakePhoton = atof(cutValue.c_str());
        }
        else if (cutName == "fakebJet")
        {
            cut.fakebJet = atof(cutValue.c_str());
        }
        else if (cutName == "topness")
        {
            cut.topness = atof(cutValue.c_str());
        }
        else if (cutName == "topness2")
        {
            cut.topness2 = atof(cutValue.c_str());
        }
        else if (cutName == "deltaPhi")
        {
            cut.deltaPhi = atof(cutValue.c_str());
        }
        else if (cutName == "HRatio")
        {
            cut.HRatio = atof(cutValue.c_str());
        }
    }

    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        event = analyseAA(branchParticle, branchTower, branchPhoton, branchElectron, branchMuon, branchJet, cut);

        if (event.status && event.hardJet.Pt() > cut.ptj/*  && event.deltaPhi > cut.deltaPhi */)
        {
            if (event.jetList.size() <= 4 /* && event.deltaPhi > 0.8 */)
            {
                //histDeltaPhi->Fill(event.deltaPhi);
                //histHRatio->Fill(max(event.higgs1.Pt(), event.higgs2.Pt()) / event.hardJet.Pt()); 
                //histThrust->Fill(event.thrust);
                histInvM->Fill(event.diHiggsInvM());
            }
            else if (event.jetList.size() > 4 && event.topness2 > cut.topness2 /* && event.deltaPhi > 1 */)
            {
                /* histNLJet->Fill(event.jetList[1].Pt() / event.jetList[0].Pt());
                histN2LJet->Fill(event.jetList[2].Pt() / event.jetList[0].Pt());
                histN3LJet->Fill(event.jetList[3].Pt() / event.jetList[0].Pt()); */
                //histTopness2->Fill(event.topness2);
                //histBJetRatio->Fill(event.b1.E() / (event.b1.E() + event.b2.E()));
                //histDeltaPhi->Fill(event.deltaPhi);
                //histHRatio->Fill(max(event.higgs1.Pt(), event.higgs2.Pt()) / event.hardJet.Pt());
                histInvM->Fill(event.diHiggsInvM());
                //histTopness->Fill(event.topness);
            }
            /* histDeltaPhi->Fill(event.deltaPhi);
            histHRatio->Fill(max(event.higgs1.Pt(), event.higgs2.Pt()) / event.hardJet.Pt()); */
        }
    } 
    //cout << argv[argc - 1] << " " << n << endl;
    TFile *fInvM = new TFile((string("../histogram/InvMass/") + string(argv[1])).c_str(), "RECREATE");
    histInvM->Write();
    fInvM->Close();
    /* TFile *fHRatio = new TFile((string("../histogram/HRatio/") + string(argv[1])).c_str(), "RECREATE");
    histHRatio->Write();
    fHRatio->Close(); */
    /* TFile *fDeltaPhi = new TFile((string("../histogram/DeltaPhi/") + string(argv[1])).c_str(), "RECREATE");
    histDeltaPhi->Write();
    fDeltaPhi->Close(); */
    /* TFile *fBJetRatio = new TFile((string("../histogram/BJetRatio/") + string(argv[1])).c_str(), "RECREATE");
    histBJetRatio->Write();
    fBJetRatio->Close(); */
    /*TFile *fNLJet = new TFile((string("../histogram/NLJet/") + string(argv[1])).c_str(), "RECREATE");
    histNLJet->Write();
    fNLJet->Close();
    TFile *fN2LJet = new TFile((string("../histogram/N2LJet/") + string(argv[1])).c_str(), "RECREATE");
    histN2LJet->Write();
    fN2LJet->Close();
    TFile *fN3LJet = new TFile((string("../histogram/N3LJet/") + string(argv[1])).c_str(), "RECREATE");
    histN3LJet->Write();
    fN3LJet->Close(); */
    /* TFile *fTopness = new TFile((string("../histogram/Topness/") + string(argv[1])).c_str(), "RECREATE");
    histTopness->Write();
    fTopness->Close();
    TFile *fTopness2 = new TFile((string("../histogram/Topness2/") + string(argv[1])).c_str(), "RECREATE");
    histTopness2->Write();
    fTopness2->Close(); */
    /*TFile *fThrust = new TFile((string("../histogram/Thrust/") + string(argv[1])).c_str(), "RECREATE");
    histThrust->Write();
    fThrust->Close();*/
    return 0;
}