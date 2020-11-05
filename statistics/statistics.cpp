#include "include/StatisticsClass.h"
#include "TFile.h"

int main()
{
    vector<int> kapLam = {0, 1, 2, 3};
    vector<double> xSection = {24.82, 16.2594, 9.276521, 8.839};
    double braRatio = 2 * 0.5809 * 0.00227;
    vector<vector<double>> sigList;
    vector<double> bkg;
    vector<vector<double>> result;

    for (int i = 0; i < kapLam.size(); i++)
    {
        char fileName[100];
        char histName[10];
        vector<double> tmp = {};
        sprintf(fileName, "/mnt/d/work/Hpair/cpp/showerAnalysis/histogram/sig_aa_%d_10_28.root", kapLam[i]);
        sprintf(histName, "kappa=%d", kapLam[i]);
        TFile f(fileName);
        TH1D *hist = (TH1D*) f.Get(histName);
        double scale = 2 * xSection[i] * braRatio * 3000 / 100000;
        tmp = histConverter(hist, scale);
        sigList.push_back(tmp);
        f.Close();
    }
    
    TFile fbkg("/mnt/d/work/Hpair/cpp/showerAnalysis/histogram/bkg_aa_10.root");
    TH1D *histBkg = (TH1D*) fbkg.Get("bkg");
    double xSecBkg = 246.3;
    bkg = histConverter(histBkg, 2 * xSecBkg * 3000 / 300000);
    fbkg.Close();
    /* double num;
    for (int i = 0; i < sigList[0].size(); i++)
    {
        num += sigList[0][i];
    }
    cout << num << endl; */

    result = sigmaCalc(sigList, sigList[1], bkg, kapLam);
    for (int i = 0; i < result.size(); i++)
    {
        cout << result[i][1] << endl;
    }
    
    
    return 0;
}