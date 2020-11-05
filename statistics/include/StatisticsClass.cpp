#include "StatisticsClass.h"

vector<double> histConverter(TH1D *hist, double scale)
{
    vector<double> tmp = {};
    for (int i = 0; i < hist->GetNbinsX(); i++)
    {
        tmp.push_back(hist->GetBinContent(i) * scale);
    }
    return tmp;
}

double logfac(int n)
{
    double result = 0;
    for (int i = 1; i <= n; i++)
    {
        result += log(i);
    }
    return result;
}

void LikelihoodFunc::init(Data data)
{
    this->signal = data.signal;
    this->bkg = data.bkg;
    this->obs = data.obs;
}

double LikelihoodFunc::likelihood(const double *sigStrength)
{
    int nBin = signal.size();
    double result = 0;
    for (int i = 0; i < nBin; i++)
    {
        double mu = sigStrength[0] * round(signal[i]) + round(bkg[i]);
        if (mu != 0)
        {
            result += (round(obs[i])) * log(mu) - logfac(round(obs[i])) - mu;
        }
    }
    return -result;
}

double maxLikelihood(Data data)
{
    LikelihoodFunc likelihood = LikelihoodFunc();
    likelihood.init(data);

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer();
    ROOT::Math::Functor f(&likelihood, &LikelihoodFunc::likelihood, 1);

    minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimizer->SetMaxIterations(10000);  // for GSL
    minimizer->SetTolerance(0.001);
    minimizer->SetPrintLevel(0);
    minimizer->SetFunction(f);

    minimizer->SetVariable(0, "sigStrength", 1.2, 0.01);
    minimizer->Minimize();
    cout << (minimizer->X())[0] << endl;
    return minimizer->MinValue();
}


vector<double> addVec(vector<double> vec1, vector<double> vec2)
{
    vector<double> sum;
    if (vec1.size() != vec2.size())
    {
        cout << "vector size does not match!" << endl;
        return {};
    }
    else
    {
        for (int i = 0; i < vec1.size(); i++)
        {
            sum.push_back(vec1[i] + vec2[i]);
        }
        return sum;
    }
}

vector<vector<double>> sigmaCalc(vector<vector<double>> sigList, vector<double> assumption, vector<double> bkg, vector<int> kapLam)
{
    vector<vector<double>> result;
    Data bestFit;
    bestFit.bkg = bkg;
    bestFit.signal = assumption;
    bestFit.obs = addVec(bkg, assumption);

    double bestFitLogL = maxLikelihood(bestFit);

    for (int i = 0; i < sigList.size(); i++)
    {
        if (sigList[i].size() != bkg.size())
        {
            cout << "signal and background histograms do not match" << endl;
            continue;
        }
        Data data;
        data.bkg = bkg;
        data.signal = sigList[i];
        data.obs = addVec(sigList[i], bkg);
        double logL = maxLikelihood(data);
        vector<double> significance = {(double) kapLam[i], sqrt(2 * abs(logL - bestFitLogL))};
        result.push_back(significance);
    }
    return result;
}