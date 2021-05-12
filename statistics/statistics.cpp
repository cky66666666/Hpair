#include "include/StatisticsClass.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "fstream"
#include "string"

TGraph* drawGraph(vector<vector<double>> result, int begin, int end)
{
    int nEntry = end - begin + 1;
    double kapLam[nEntry], significance[nEntry];

    for (int i = begin; i <= end; i++)
    {
        kapLam[i - begin] = result[i][0];
        significance[i - begin] = result[i][1];
    }

    TGraph *graph = new TGraph(nEntry, kapLam, significance);
    return graph;
}

vector<vector<double>> scaleSignal(vector<vector<double>> sigList, int sm)
{
    vector<double> eventNum;
    vector<vector<double>> sigListScaled;

    for (int i = 0; i < sigList.size(); i++)
    {
        double num = 0;
        for (int j = 0; j < sigList[i].size(); j++)
        {
            num += sigList[i][j];
        }
        eventNum.push_back(num);
    }
    
    for (int i = 0; i < sigList.size(); i++)
    {
        vector<double> histScaled;
        double scale = eventNum[sm] / eventNum[i];
        for (int j = 0; j < sigList[i].size(); j++)
        {
            histScaled.push_back(sigList[i][j] * scale);
        }
        sigListScaled.push_back(histScaled);
    }

    return sigListScaled;
}

vector<double> bkgCombiner(int nBin)
{
    TFile *fbkg_tth = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_tth_100_100.root");
    TFile *fbkg_tthj = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_tthj_100_100.root");
    TFile *fbkg_bbaa = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_bbaa_100_100.root");
    TFile *fbkg_bbaj = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_bbaj_100_100.root");
    TFile *fbkg_bbjj = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_bbjj_100_100.root");
    TFile *fbkg_bjaa = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_bjaa_100_100.root");
    TFile *fbkg_jjaa = new TFile("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_jjaa_100_100.root");

    const double lumi = 30000;
    const double braRatioAA = 2 * 0.5809 * 0.00227;
    const double fake = 0.0005;

    vector<TFile*> bkgA = {fbkg_tth, fbkg_tthj, fbkg_bbaa, fbkg_bbaj, fbkg_bbjj, fbkg_bjaa, fbkg_jjaa};
    vector<string> bkgNameA = {"bkg_tth_InvMass", "bkg_tthj_InvMass", "bkg_bbaa_InvMass", "bkg_bbaj_InvMass", "bkg_bbjj_InvMass", "bkg_bjaa_InvMass", "bkg_jjaa_InvMass"};
    vector<double> bkgScaleA = {4.396 * lumi / 1000000, 30.68 * lumi / 1000000, 45.13 * lumi / 1000000, 2 * fake * 252000 * lumi / 1000000, 
    fake * fake * 70000000 * lumi / 1000000, 2985.6542 * 0.025 * lumi / 1000000, 14477.378 * 0.025 * 0.025 * lumi / 1000000};

    vector<double> bkg = {};
    vector<double> bkgN = {};
    TH1D *bkgHist = new TH1D("bkg", "bkg", 50, 250, 1750);

    for (int i = 0; i < bkgA.size(); i++)
    {
        char histName[bkgNameA[i].length()];
        strcpy(histName, bkgNameA[i].c_str());
        TH1D *hist = (TH1D*) bkgA[i]->Get(histName);
        bkgHist->Add(hist, bkgScaleA[i]);
        bkgN.push_back(hist->GetEntries() * bkgScaleA[i]);
    }
    for (int i = 1; i <= nBin; i++)
    {
        bkg.push_back(bkgHist->GetBinContent(i));
    }

    delete bkgHist;
    return bkg;
}

void output(vector<double> bkg, vector<double> signal)
{
    fstream file("cut.txt", ios::app);
    string content;
    double bkgN = 0;

    for (int i = 0; i < bkg.size(); i++)
    {
        char tmp[20];
        sprintf(tmp, "%.2lf", bkg[i]);
        content.append(string(tmp) + " ");
        bkgN += bkg[i];
    }
    
    char tmp[20];
    sprintf(tmp, "%.2lf", signal[3] / sqrt(bkgN));
    content.append(string(tmp));
    file.seekg(0, ios::end);
    file << content;
    file << endl;
}

int main()
{
    vector<float> kapLam = {-2, -1, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5};
    vector<string> fileIndex = {"-2", "-1", "0", "0_2", "0_4", "0_6", "0_8", "1", "1_2", "1_4", "1_6", "1_8", "2", "3", "4", "5"};
    vector<string> histIndex = {"-2", "-1", "0", "0.2", "0.4", "0.6", "0.8", "1", "1.2", "1.4", "1.6", "1.8", "2", "3", "4", "5"};
    vector<double> xSection = {747.2, 447.4, 277.1, 245.9, 217.7, 191.4, 168.4, 147.9, 130.5, 115.7, 103.7, 94.78, 88.4, 98.73, 179.3, 329.8};
    const double braRatioAA = 2 * 0.5809 * 0.00227;
    const double braRatioTT = 2 * 0.5809 * 0.06256;
    const double lumi = 30000;
    const double mcNum = 100000;
    vector<vector<double>> sigListAA, sigListTT;
    vector<double> bkgAA, bkgTT, sigAA;
    vector<vector<double>> resultAA, resultTT, resultComb;
    
    for (int i = 0; i < fileIndex.size(); i++)
    {
        char fileNameAA[200], fileNameTT[200];
        char histName[50];
        vector<double> tmp = {};
        sprintf(fileNameAA, "/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/sig_aa_%s_10_100.root", fileIndex[i].c_str());
        //cout << fileNameAA << endl;
        sprintf(histName, "kappa=%s_InvMass", histIndex[i].c_str());
        //cout << histName << endl;
        TFile fA(fileNameAA);
        TH1D *histA = (TH1D*) fA.Get(histName);
        //cout << hist->GetEntries() << endl;
        double scaleA = xSection[i] * braRatioAA * lumi / mcNum;
        double scaleT = xSection[i] * braRatioTT * lumi / mcNum;
        tmp = histConverter(histA, scaleA, 50);
        sigListAA.push_back(tmp);
        sigAA.push_back(histA->GetEntries() * scaleA);
        fA.Close();
        //fT.Close();
    }

    double num;
    
    //sigListAA = scaleSignal(sigListAA, 3);
    bkgAA = bkgCombiner(50);
    resultAA = sigmaCalc(sigListAA, sigListAA[7], bkgAA, kapLam);
    /* cout << resultAA[2][1] << endl;
    for (int i = 0; i < resultAA.size(); i++)
    {
        cout << resultAA[i][1] << endl;
    } */
    
    //cout << resultAA[2][0] << endl;
    //output(bkgAA, sigAA);
    //resultTT = sigmaCalc(sigListTT, sigListTT[2], bkgTT, kapLam);

    /* for (int i = 0; i < resultAA.size(); i++)
    {
        resultComb.push_back({(double) kapLam[i], sqrt(resultAA[i][1] * resultAA[i][1] + resultTT[i][1] * resultTT[i][1])});
    } */
    
    
    TGraph *graphAA = drawGraph(resultAA, 2, 12);
    TFile fAA("../graph/AA_100_120_All.root", "RECREATE");
    graphAA->Write();
    fAA.Close();

    /* TGraph *graphTT = drawGraph(resultTT);
    TFile fTT("../graph/significanceTT.root", "RECREATE");
    graphTT->Write();
    fTT.Close(); */

    /* TGraph *graphComb = drawGraph(resultComb);
    TFile fComb("../graph/significance.root", "RECREATE");
    graphComb->Write();
    fComb.Close(); */
    
    return 0;
}