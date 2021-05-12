#ifndef _StatisticsClass_
#define _StatisticsClass_

#include "vector"
#include "TH1D.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

using namespace std;

struct Data
{
    vector<double> signal, bkg, obs; 
};

vector<double> histConverter(TH1D*, double, int);
double logfac(int);
double maxLikelihood(Data);
vector<double> addVec(vector<double>, vector<double>);
vector<vector<double>> sigmaCalc(vector<vector<double>>, vector<double>, vector<double>, vector<float>);

class LikelihoodFunc
{
private:
    vector<double> signal, bkg, obs;

public:

    LikelihoodFunc()
    {

    }
    ~LikelihoodFunc()
    {

    }

    void init(Data);
    double likelihood(const double*);
};



#endif