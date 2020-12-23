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

int main()
{
    vector<int> kapLam = {-1, 0, 1, 2, 3, 4, 5};
    vector<double> xSection = {44.01856, 24.82, 16.2594, 9.276521, 8.839, 17.25444, 32.60771};
    const double braRatioAA = 2 * 0.5809 * 0.00227;
    const double braRatioTT = 2 * 0.5809 * 0.06256;
    vector<vector<double>> sigListAA, sigListTT;
    vector<double> bkgAA, bkgTT;
    vector<vector<double>> resultAA, resultTT, resultComb;

    for (int i = 0; i < kapLam.size(); i++)
    {
        char fileNameAA[200], fileNameTT[200];
        char histName[10];
        vector<double> tmp = {};
        sprintf(fileNameAA, "/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/sig_aa_%d_10_28.root", kapLam[i]);
        sprintf(fileNameTT, "/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/sig_tt_%d_10_28.root", kapLam[i]);
        sprintf(histName, "kappa=%d", kapLam[i]);
        TFile fA(fileNameAA)/* , fT(fileNameTT) */;
        TH1D *histA = (TH1D*) fA.Get(histName);
        //TH1D *histT = (TH1D*) fT.Get(histName);
        //cout << hist->GetEntries() << endl;
        double scaleA = xSection[i] * braRatioAA * 3000 / 100000;
        double scaleT = xSection[i] * braRatioTT * 3000 / 100000;
        tmp = histConverter(histA, scaleA);
        sigListAA.push_back(tmp);
        /* tmp = histConverter(histT, scaleT);
        sigListTT.push_back(tmp); */
        fA.Close();
        //fT.Close();
    }

    
    TFile fbkgAA("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_28.root");
    //TFile fbkgTT("/home/E/chaikangyu/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_tt_28.root");
    TH1D *histBkgAA = (TH1D*) fbkgAA.Get("bkg");
    //TH1D *histBkgTT = (TH1D*) fbkgTT.Get("bkg");
    bkgAA = histConverter(histBkgAA, 1);
    for (int i = 0; i < sigListAA[0].size(); i++)
    {
        cout << sigListAA[0][i] << endl;
    }
    
    //bkgTT = histConverter(histBkgTT, 1);
    fbkgAA.Close();
    //fbkgTT.Close();
    double num;
    /* for (int i = 0; i < bkg.size(); i++)
    {
        num += bkg[i];
    }
    cout << num << endl; */

    resultAA = sigmaCalc(sigListAA, sigListAA[2], bkgAA, kapLam);
    //resultTT = sigmaCalc(sigListTT, sigListTT[2], bkgTT, kapLam);

    /* for (int i = 0; i < resultAA.size(); i++)
    {
        resultComb.push_back({(double) kapLam[i], sqrt(resultAA[i][1] * resultAA[i][1] + resultTT[i][1] * resultTT[i][1])});
    } */
    
    
    TGraph *graphAA = drawGraph(resultAA);
    TFile fAA("../graph/significanceAA.root", "RECREATE");
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