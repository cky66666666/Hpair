#include "iostream"
#include "string.h"
#include "vector"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootClasses.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "algorithm"

using namespace std;


vector<TRootLHEFParticle*> findParticles(TClonesArray *branchParticles){
    vector<TRootLHEFParticle*> particleList;
    TRootLHEFParticle *particle;
    for (int i = 0; i < branchParticles->GetEntries(); i++)
    {
        particle = (TRootLHEFParticle*) branchParticles->At(i);
        if (particle->Status == 1)
        {
            particleList.push_back(particle);
        }
    }
    return particleList;
}

TH1D* drawHist(ExRootTreeReader *treeReader, char treeName[], int ptj){
    TClonesArray *branchParticles = treeReader->UseBranch("Particle");
    TRootLHEFParticle *particle;
    vector<TLorentzVector> higgs, photon, bQuark;
    TLorentzVector tmp;
    int maxPtPoint, minPtPoint;
    /* char histName[10];
    sprintf(histName, "ptj=%d", ptj); */

    TH1D *histInvMass = new TH1D(treeName, treeName, 50, 250, 1000);
    //TH1D *histRH = new TH1D(histName, histName, 50, 0, 5);
    for (int iEvent = 0; iEvent < treeReader->GetEntries(); iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        for (int iParticle = 0; iParticle < branchParticles->GetEntries(); iParticle++)
        {
            particle = (TRootLHEFParticle*) branchParticles->At(iParticle);
            tmp.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
            if (particle->PID == 25)
            {
                higgs.push_back(tmp);
            }
            else if (particle->PID == 22)
            {
                photon.push_back(tmp);
            }
            else if (abs(particle->PID) == 5)
            {
                bQuark.push_back(tmp);
            }
            
        }
        //int maxPoint = maxPt(higgs);
        //cout << photon.size() << endl;
        if (higgs.size() != 2 || photon.size() != 2 || bQuark.size() != 2)
        {
            higgs.clear();
            photon.clear();
            bQuark.clear();
            continue;
        }
        if ((higgs[0] + higgs[1]).Pt() > ptj && photon[0].Pt() > 30 && photon[1].Pt() > 30 && bQuark[0].Pt() > 30 && bQuark[1].Pt() > 30
            && bQuark[0].DeltaR(bQuark[1]) < 2 && bQuark[0].DeltaR(bQuark[1]) > 0.4 && photon[0].DeltaR(photon[1]) > 0.4 && photon[0].DeltaR(photon[1]) < 2
            && bQuark[0].DeltaR(photon[0]) > 0.4 && bQuark[0].DeltaR(photon[1]) > 0.4 && bQuark[1].DeltaR(photon[0]) > 0.4 && bQuark[1].DeltaR(photon[1]) > 0.4)
        {
           histInvMass->Fill((higgs[0] + higgs[1]).M());
        }

        //histRH->Fill(higgs[0].DeltaR(higgs[1]));
        higgs.clear();
        photon.clear();
        bQuark.clear();
    }
    //histInvMass->Scale(10000 / (histInvMass->GetEntries()), "nosw2");
    //cout << n << endl;
    cout << histInvMass->GetEntries() << endl;
    return histInvMass;
}

vector<TLorentzVector> findHiggs(vector<TLorentzVector> higgsCand)
{
    TLorentzVector h1, h2;
    vector<TLorentzVector> tmp;
    h1.SetPxPyPzE(0, 0, 0, 0);
    h2.SetPxPyPzE(0, 0, 0, 0);

    if (higgsCand.size() < 4) return {h1, h2};

    for (int i = 0; i < higgsCand.size(); i++)
    {
        for (int j = i + 1; j < higgsCand.size(); j++)
        {
            tmp = higgsCand;
            tmp.erase(tmp.begin() + j);
            tmp.erase(tmp.begin() + i);
            if (abs((higgsCand[i] + higgsCand[j]).M() - 125) < 20 && abs((tmp[0] + tmp[1]).M() - 125) < 20)
            {
                h1 = higgsCand[i] + higgsCand[j];
                h2 = tmp[0] + tmp[1];
                i = 100;
                j = 100;
            }
            
        }
        
    }
    
    return {h1, h2};
}

TH1D* analyzeBkg(ExRootTreeReader *treeReader)
{
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    int nEvent = treeReader->GetEntries(), n = 0;
    TRootLHEFParticle *particle;
    vector<TLorentzVector> higgsCand, higgs;
    TLorentzVector pb;
    TH1D *histBkg = new TH1D("bkg", "bkg", 50, 0, 1000);

    for (int iEvent = 0; iEvent < nEvent; iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        for (int i = 0; i < branchParticle->GetEntriesFast(); i++)
        {
            particle = (TRootLHEFParticle*) branchParticle->At(i);
            if (abs(particle->PID) == 5)
            {
                pb.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
                higgsCand.push_back(pb);
            }
        }
        
        higgs = findHiggs(higgsCand);
        /* if (higgs[0].Pt() > higgsPt && higgs[1].Pt() > higgsPt && (higgs[0] + higgs[1]).Pt() > jetPt)
        {
            histBkg->Fill(higgs[0].Pt());
            histBkg->Fill(higgs[1].Pt());
            n += 1;
        } */
        higgsCand.clear();
    }
    
    cout << n << endl;
    return histBkg;
}

int main(int argc, char *argv[]){
    vector<string> treeName_s = {"kappa=1", "kappa=2", "kappa=3"};
    vector<string> inputFile_s = {"../event/sig_aa_1_10_28.root", "../event/sig_aa_2_10_28.root", "../event/sig_aa_3_10_28.root"};
    vector<int> nEvent;
    THStack *stackInvMass = new THStack("InvMass", "InvMass");
    for (int i = 0; i < treeName_s.size(); i++)
    {
        char treeName[treeName_s[i].length()], inputFile[inputFile_s[i].length()];
        strcpy(treeName, treeName_s[i].c_str());
        strcpy(inputFile, inputFile_s[i].c_str());

        TChain *chain = new TChain(treeName);
        chain->Add(inputFile);
        ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
        stackInvMass->Add(drawHist(treeReader, treeName, 200));
        //stackInvMass->Add(drawHist(treeReader, treeName, 150));
        //stackInvMass->Add(drawHist(treeReader, treeName, 200));
        //nEvent.push_back(drawHist(treeReader, nTree[i]));
        
        
        delete chain;
        delete treeReader;
    } 

    /* char bkg[] = "bkg4b";
    TChain *chain = new TChain(bkg);
    chain->Add("../bkg4b.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    stackInvMass->Add(analyzeBkg(treeReader)); */
    //int nBkg = analyzeBkg(treeReader);
    /* for (int i = 0; i < nEvent.size(); i++)
    {
        cout << nEvent[i] / nBkg << endl;
    } */
    
    TFile f("../hist.root", "RECREATE");
    stackInvMass->Write();
    f.Close();
    return 0;
}