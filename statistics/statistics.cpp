#include "include/StatisticsClass.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"

TGraph* drawGraph(vector<vector<double>> result)
{
    int nEntry = result.size();
    double kapLam[nEntry], significance[nEntry];
    for (int i = 0; i < nEntry; i++)
    {
        kapLam[i] = result[i][0];
        significance[i] = result[i][1];
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

int main()
{
    vector<int> kapLam = {-2, -1, 0, 1, 2, 3, 4, 5};
    vector<double> xSection = {747.2, 447.4, 277.1, 147.9, 88.4, 98.73, 179.3, 329.8};
    const double braRatioAA = 2 * 0.5809 * 0.00227;
    const double braRatioTT = 2 * 0.5809 * 0.06256;
    const double lumi = 30000;
    const double mcNum = 100000;
    vector<vector<double>> sigListAA, sigListTT;
    vector<double> bkgAA, bkgTT;
    vector<vector<double>> resultAA, resultTT, resultComb;

    for (int i = 0; i < kapLam.size(); i++)
    {
        char fileNameAA[200], fileNameTT[200];
        char histName[10];
        vector<double> tmp = {};
        sprintf(fileNameAA, "/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/sig_aa_%d_10_100.root", kapLam[i]);
        sprintf(fileNameTT, "/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/sig_tt_%d_10_100.root", kapLam[i]);
        sprintf(histName, "kappa=%d", kapLam[i]);
        TFile fA(fileNameAA)/* , fT(fileNameTT) */;
        TH1D *histA = (TH1D*) fA.Get(histName);
        //TH1D *histT = (TH1D*) fT.Get(histName);
        //cout << hist->GetEntries() << endl;
        double scaleA = xSection[i] * braRatioAA * lumi / mcNum;
        double scaleT = xSection[i] * braRatioTT * lumi / mcNum;
        tmp = histConverter(histA, scaleA, 30);
        sigListAA.push_back(tmp);
        /* tmp = histConverter(histT, scaleT);
        sigListTT.push_back(tmp); */
        fA.Close();
        //fT.Close();
    }

    
    TFile fbkgAA("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_100.root");
    //TFile fbkgTT("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_tt_28.root");
    TH1D *histBkgAA = (TH1D*) fbkgAA.Get("bkg");
    //TH1D *histBkgTT = (TH1D*) fbkgTT.Get("bkg");
    bkgAA = histConverter(histBkgAA, 1, 30);
    /* for (int i = 0; i < sigListAA[0].size(); i++)
    {
        cout << sigListAA[0][i] << endl;
    } */
    
    //bkgTT = histConverter(histBkgTT, 1);
    fbkgAA.Close();
    //fbkgTT.Close();
    double num;
    
    sigListAA = scaleSignal(sigListAA, 3);
    resultAA = sigmaCalc(sigListAA, sigListAA[3], bkgAA, kapLam);
    //resultTT = sigmaCalc(sigListTT, sigListTT[2], bkgTT, kapLam);

    /* for (int i = 0; i < resultAA.size(); i++)
    {
        resultComb.push_back({(double) kapLam[i], sqrt(resultAA[i][1] * resultAA[i][1] + resultTT[i][1] * resultTT[i][1])});
    } */
    
    
    TGraph *graphAA = drawGraph(resultAA);
    TFile fAA("../graph/AA_100_all_scaled.root", "RECREATE");
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