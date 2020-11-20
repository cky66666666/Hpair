#include "include/StatisticsClass.h"
#include "TFile.h"
#include "TGraph.h"

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
    double braRatio = 2 * 0.5809 * 0.00227;
    vector<vector<double>> sigList;
    vector<double> bkg;
    vector<vector<double>> result;

    for (int i = 0; i < kapLam.size(); i++)
    {
        char fileName[100];
        char histName[10];
        vector<double> tmp = {};
        sprintf(fileName, "/mnt/d/work/Hpair/cpp/showerAnalysis/histogram/InvMass/sig_aa_%d_10_28.root", kapLam[i]);
        sprintf(histName, "kappa=%d", kapLam[i]);
        TFile f(fileName);
        TH1D *hist = (TH1D*) f.Get(histName);
        double scale = xSection[i] * braRatio * 3000 / 100000;
        tmp = histConverter(hist, scale);
        sigList.push_back(tmp);
        f.Close();
    }
    
    TFile fbkg("/mnt/d/work/Hpair/cpp/showerAnalysis/histogram/InvMass/bkg_aa_28.root");
    TH1D *histBkg = (TH1D*) fbkg.Get("bkg");
    bkg = histConverter(histBkg, 1);
    fbkg.Close();
    double num;
    /* for (int i = 0; i < bkg.size(); i++)
    {
        num += bkg[i];
    }
    cout << num << endl; */

    result = sigmaCalc(sigList, sigList[2], bkg, kapLam);
    for (int i = 0; i < result.size(); i++)
    {
        cout << result[i][1] << endl;
    }
    TGraph *graph = drawGraph(result);
    TFile f("../graph/significance.root", "RECREATE");
    graph->Write();
    f.Close();
    
    return 0;
}