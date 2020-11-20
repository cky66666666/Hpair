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

TH1D* drawHist(ExRootTreeReader *treeReader, char treeName[], int ptj, double scale)
{
    TClonesArray *branchParticles = treeReader->UseBranch("Particle");
    TRootLHEFParticle *particle;
    vector<TLorentzVector> higgs, photon, bQuark;
    vector<double> jetPt;
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
            if (particle->PID == 22 && particle->PT > 30)
            {
                photon.push_back(tmp);
            }
            else if (abs(particle->PID) == 5 && particle->PT > 30)
            {
                bQuark.push_back(tmp);
            }
            else if (abs(particle->PID) <= 4 || abs(particle->PID) == 21)
            {
                jetPt.push_back(tmp.Pt());
            }
        }
        if (photon.size() != 2 || bQuark.size() != 2)
        {
            higgs.clear();
            photon.clear();
            bQuark.clear();
            jetPt.clear();
            continue;
        }
        else
        {
            higgs.push_back(photon[0] + photon[1]);
            higgs.push_back(bQuark[0] + bQuark[1]);
        }

        int maxPoint = max_element(jetPt.begin(), jetPt.end()) - jetPt.begin();

        if (jetPt[maxPoint] > ptj && photon[0].Pt() > 30 && photon[1].Pt() > 30 && bQuark[0].Pt() > 30 && bQuark[1].Pt() > 30
            && bQuark[0].DeltaR(bQuark[1]) < 2 && bQuark[0].DeltaR(bQuark[1]) > 0.4 && photon[0].DeltaR(photon[1]) > 0.4 && photon[0].DeltaR(photon[1]) < 2
            && bQuark[0].DeltaR(photon[0]) > 0.4 && bQuark[0].DeltaR(photon[1]) > 0.4 && bQuark[1].DeltaR(photon[0]) > 0.4 && bQuark[1].DeltaR(photon[1]) > 0.4
            && abs(higgs[0].M() - 125) < 3 && abs(higgs[1].M() - 125) < 20)
        {
           histInvMass->Fill((higgs[0] + higgs[1]).M());
        }

        //histRH->Fill(higgs[0].DeltaR(higgs[1]));
        higgs.clear();
        photon.clear();
        bQuark.clear();
        jetPt.clear();
    }
    //histInvMass->Scale(10000 / (histInvMass->GetEntries()), "nosw2");
    //cout << n << endl;
    cout << histInvMass->GetEntries() << endl;
    histInvMass->Scale(scale, "nosw2");
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
    vector<string> treeName_s = {"kappa=3", "bkg_3jh", "bkg_tth", "bkg_zh"};
    vector<string> histName_s = {"kappa=3", "bkg_3jh", "bkg_tth", "bkg_zh"};
    vector<string> inputFile_s = {"../event/sig_aa_3_10_28.root", "../event/bkg_aa_3jh_20_28.root", "../event/bkg_aa_tth_10_28.root"
                                    , "../event/bkg_aa_zh_10_28.root"};
    const double braRatio = 2 * 0.5809 * 0.00227;
    const double lumi = 3000;
    vector<double> xSection = {8.839 * braRatio * lumi / 1000, 0.06625 * lumi / 20000, 2.319 * lumi / 10000, 0.05442 * lumi / 10000};
    vector<int> nEvent;
    THStack *stackInvMass = new THStack("InvMass", "InvMass");
    for (int i = 0; i < treeName_s.size(); i++)
    {
        char treeName[treeName_s[i].length()], inputFile[inputFile_s[i].length()], histName[histName_s[i].length()];
        strcpy(treeName, treeName_s[i].c_str());
        strcpy(inputFile, inputFile_s[i].c_str());
        strcpy(histName, histName_s[i].c_str());

        TChain *chain = new TChain(treeName);
        chain->Add(inputFile);
        ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
        stackInvMass->Add(drawHist(treeReader, histName, 300, xSection[i]));
        /* stackInvMass->Add(drawHist(treeReader, treeName, 200));
        stackInvMass->Add(drawHist(treeReader, treeName, 300)); */
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