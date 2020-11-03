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
#include "Tagger/BAChannel.h"
#include "Tagger/BTChannel.h"

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
    if (abs(event.softHiggs.M() - 125) > 20 || event.hardJet.Pt() < 150 || event.hardHiggs.Pt() < 80 || event.softHiggs.Pt() < 80)
    {
        status = false;
    }
    return status;
}

double analyseBB(TClonesArray *branchJet, TClonesArray *branchParticle, TClonesArray *branchTower, int type)
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

        boostedHiggs->init(finalState, parton, {}, 80, 110, 1.5, 0.3);
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

double analyseAA(TClonesArray *branchParticle, TClonesArray *branchTower, TClonesArray *branchPhoton, TClonesArray *branchElectron, TClonesArray *branchMuon, TClonesArray *branchJet)
{
    BAChannel *BAEvent = new BAChannel();
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
    

    BAEvent->init(finalState, parton, photon, electron, muon, delphesJet);
    BAEvent->process();

    if (BAEvent->status && (BAEvent->hardJet).Pt() > 200 /* && (BAEvent->higgsFromA).Pt() > 80 && (BAEvent->higgsFromB).Pt() > 80 */)
    {
        inv = (BAEvent->higgsFromA + BAEvent->higgsFromB).M();
        BAEvent->finish();
        delete BAEvent;
        return inv;
    }
    else
    {
        BAEvent->finish();
        delete BAEvent;
        return 0;
    }
}

double analyzeTT(TClonesArray *branchJet, TClonesArray *branchElectron, TClonesArray *branchMuon, TClonesArray *branchParticle, TClonesArray *branchMET)
{
    vector<Jet*> jet;
    vector<Electron*> electron;
    vector<Muon*> muon;
    vector<GenParticle*> parton;
    TLorentzVector missP;

    for (int i = 0; i < branchJet->GetEntries(); i++)
    {
        jet.push_back((Jet*) branchJet->At(i));
    }
    for (int i = 0; i < branchElectron->GetEntries(); i++)
    {
        electron.push_back((Electron*) branchElectron->At(i));
    }
    for (int i = 0; i < branchMuon->GetEntries(); i++)
    {
        muon.push_back((Muon*) branchMuon->At(i));
    }
    MissingET *met = (MissingET*) branchMET->At(0);
    missP.SetPtEtaPhiM(met->MET, met->Eta, met->Phi, 0);
    
    BTChannel BTEvent = BTChannel();
    BTEvent.init(jet, electron, muon, missP);
    BTEvent.process();

    if (BTEvent.status)
    {
        ROOT::Math::Minimizer *mt2 = ROOT::Math::Factory::CreateMinimizer();
        ROOT::Math::Functor f(&BTEvent, &BTChannel::mTMax, 2);
        TRandom2 r(5);
        mt2->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
        mt2->SetMaxIterations(10000);  // for GSL
        mt2->SetTolerance(0.001);
        mt2->SetPrintLevel(0);
        mt2->SetFunction(f);

        mt2->SetVariable(0, "x1", 10, 0.01);
        mt2->SetVariable(1, "x2", 10, 0.01);

        mt2->Minimize();
        if (mt2->MinValue() > 100)
        {
            return BTEvent.dihiggsInvM;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
    
}

int main(int argc, char *argv[])
{
    // Usage: ./eventAnalysis /path/to/hist.root /path/to/inputfile histogram name

    //gSystem->Load("/mnt/d/work/Hpair/Delphes/libDelphes");
    
    TH1D *hist = new TH1D(argv[3], argv[3], 50, 250, 1000);
    TFile *f = new TFile(argv[1], "RECREATE");

    TChain *chain = new TChain("Delphes");
    chain->Add(argv[2]);
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");

    //cout << treeReader->GetEntries() << endl;
    int nEvent = treeReader->GetEntries();
    int n = 0;
    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        //cout << iEvent << endl;
        //double inv = analyseBB(branchJet, branchParticle, branchTower, 1);
        double inv = analyseAA(branchParticle, branchTower, branchPhoton, branchElectron, branchMuon, branchJet);
        //double inv = analyzeTT(branchJet, branchElectron, branchMuon, branchParticle, branchMET);
        if (inv > 0)
        {
            hist->Fill(inv);
            //cout << inv << endl;
            n += 1;
        }
        treeReader->Clear();
    }
    cout << n << endl;
    hist->Write();    
    f->Close();

    return 0;
}