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

#define MH 125
#define PTJ 100
#define NBIN 50
#define BEGIN 250
#define END 1000

using namespace std;

struct MyEvent
{
    TLorentzVector bJet1, bJet2, photon1, photon2, hardJet;
    MyEvent()
    {
        bJet1.SetPxPyPzE(0, 0, 0, 0);
        bJet2.SetPxPyPzE(0, 0, 0, 0);
        photon1.SetPxPyPzE(0, 0, 0, 0);
        photon2.SetPxPyPzE(0, 0, 0, 0);
        hardJet.SetPxPyPzE(0, 0, 0, 0);
    }
    ~MyEvent()
    {

    }
};

TLorentzVector findMaxPt(vector<TLorentzVector> jet)
{
    vector<double> pt = {};
    for (int i = 0; i < jet.size(); i++)
    {
        pt.push_back(jet[i].Pt());
    }
    int posMax = max_element(pt.begin(), pt.end()) - pt.begin();
    return jet[posMax];
}

vector<int> findHiggsB(vector<TLorentzVector> bJet)
{
    vector<int> index = {0, 1};
    double deltaInvMass = 1000;
    for (int i = 0; i < bJet.size(); i++)
    {
        for (int j = i + 1; j < bJet.size(); j++)
        {
            if (abs((bJet[i] + bJet[j]).M() - MH) < deltaInvMass)
            {
                deltaInvMass = abs((bJet[i] + bJet[j]).M() - MH);
                index[0] = i;
                index[1] = j;
            }
        }
        
    }
    return index;
}

MyEvent eventCons(TClonesArray *branchParticle)
{
    vector<TLorentzVector> bJet, photon, jet;
    TRootLHEFParticle *particle;
    TLorentzVector tmp;
    MyEvent event;

    for (int iParticle = 0; iParticle < branchParticle->GetEntries(); iParticle++)
    {
        particle = (TRootLHEFParticle*) branchParticle->At(iParticle);
        tmp.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
        if (particle->PID == 22 && particle->PT >= 30)
        {
            photon.push_back(tmp);
        }
        else if (abs(particle->PID) == 5 && particle->PT >= 30)
        {
            bJet.push_back(tmp);
        }
        else if (abs(particle->PID) <= 4 || particle->PID == 21)
        {
            jet.push_back(tmp);
        }
    }

    if (photon.size() != 2 || bJet.size() < 2 || jet.size() == 0)
    {
        return event;
    }
    else
    {
        if (bJet.size() == 2)
        {
            event.bJet1 = bJet[0];
            event.bJet2 = bJet[1];
            event.photon1 = photon[0];
            event.photon2 = photon[1];
            event.hardJet = findMaxPt(jet);
        }
        else
        {
            vector<int> index = findHiggsB(jet);
            event.photon1 = photon[0];
            event.photon2 = photon[1];
            event.bJet1 = bJet[index[0]];
            event.bJet2 = bJet[index[1]];
            bJet.erase(bJet.begin() + index[0]);
            bJet.erase(bJet.begin() + index[1]);

            jet.reserve(bJet.size() + jet.size());
            jet.insert(jet.end(), bJet.begin(), bJet.end());
            event.hardJet = findMaxPt(jet);
        }
        return event;
    }
}


TH1D* drawHist(ExRootTreeReader *treeReader, char histName[], double scale)
{
    TClonesArray *branchParticles = treeReader->UseBranch("Particle");
    MyEvent event;
    /* char histName[10];
    sprintf(histName, "ptj=%d", ptj); */

    TH1D *histInvMass = new TH1D(histName, histName, NBIN, BEGIN, END);
    //TH1D *histRH = new TH1D(histName, histName, 50, 0, 5);
    //cout << treeReader->GetEntries() << endl;
    TRootLHEFParticle *particle;
    for (int iEvent = 0; iEvent < treeReader->GetEntries(); iEvent++)
    {
        treeReader->ReadEntry(iEvent);
        //particle = (TRootLHEFParticle*) branchParticles->At(1);
        //cout << particle->PT << endl;
        event = eventCons(branchParticles);
        if (abs((event.photon1 + event.photon2).M() - MH) < 3 && abs((event.bJet1 + event.bJet2).M() - MH) < 25 && event.hardJet.Pt() > PTJ)
        {
            histInvMass->Fill((event.photon1 + event.photon2 + event.bJet1 + event.bJet2).M());
            //cout << (event.photon1 + event.photon2 + event.bJet1 + event.bJet2).M() << endl;
        }
    }
    histInvMass->Scale(scale, "nosw2");
    return histInvMass;
}

TH1D* bkgAnalyzer(vector<string> bkgTreeName_s, vector<string> bkgInputFile_s, vector<double> xSectionBkg)
{
    TH1D *histInvMass = new TH1D("bkg", "bkg", NBIN, BEGIN, END);
    
    for (int i = 0; i < bkgTreeName_s.size(); i++)
    {
        char inputFile[bkgInputFile_s[i].length()], treeName[bkgTreeName_s[i].length()];
        char histName[5];
        strcpy(inputFile, bkgInputFile_s[i].c_str());
        strcpy(treeName, bkgTreeName_s[i].c_str());
        strcpy(histName, ((string) "bkg").c_str());
        
        TChain *chain = new TChain(treeName);
        chain->Add(inputFile);
        ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

        histInvMass->Add(drawHist(treeReader, histName, xSectionBkg[i]));
    }

    return histInvMass;
}


int main(int argc, char *argv[])
{
    vector<string> sigTreeName_s = {"kappa=-1", "kappa=1", "kappa=3"};
    vector<string> sigHistName_s = {"kappa=-1", "kappa=1", "kappa=3"};
    vector<string> sigInputFile_s = {"../events/sig_aa_-1_10_28.root", "../events/sig_aa_1_10_28.root", "../events/sig_aa_3_10_28.root"};

    vector<string> bkgInputFile_s = {"../events/bkg_aa_3jh_20_28.root", "../events/bkg_aa_bbaa_10_28.root", 
                                    "../events/bkg_aa_tth_10_28.root", "../events/bkg_aa_zh_10_28.root"};
    vector<string> bkgTreeName_s = {"bkg_3jh", "bkg_bbaa", "bkg_tth", "bkg_zh"};

    const double braRatio = 2 * 0.5809 * 0.00227;
    const double lumi = 3000;
    vector<double> xSectionSignal = {44.01856 * braRatio * lumi / 10000, 16.2594 * lumi / 10000, 8.839 * lumi / 10000};
    vector<double> xSectionBkg = {0.06625 * lumi / 200000, 0.44785340233 * lumi / 100000, 2.319 * lumi / 100000, 0.05442 * lumi / 100000};

    vector<int> nEvent;
    THStack *stackInvMass = new THStack("InvMass", "InvMass");
    for (int i = 0; i < sigTreeName_s.size(); i++)
    {
        char treeName[sigTreeName_s[i].length()], inputFile[sigInputFile_s[i].length()], histName[sigHistName_s[i].length()];
        strcpy(treeName, sigTreeName_s[i].c_str());
        strcpy(inputFile, sigInputFile_s[i].c_str());
        strcpy(histName, sigHistName_s[i].c_str());

        TChain *chain = new TChain(treeName);
        chain->Add(inputFile);
        ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
        stackInvMass->Add(drawHist(treeReader, histName, xSectionSignal[i]));
        /* stackInvMass->Add(drawHist(treeReader, treeName, 200));
        stackInvMass->Add(drawHist(treeReader, treeName, 300)); */
        //nEvent.push_back(drawHist(treeReader, nTree[i]));
        delete chain;
        delete treeReader;
    } 

    stackInvMass->Add(bkgAnalyzer(bkgTreeName_s, bkgInputFile_s, xSectionBkg));
    
    TFile f("../hist/hist.root", "RECREATE");
    stackInvMass->Write();
    f.Close();
    return 0;
}