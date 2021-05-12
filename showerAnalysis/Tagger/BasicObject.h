#ifndef _BasicObject_
#define _BasicObject_

#include "vector"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

using namespace std;


template<typename TYPE>
vector<TYPE> combineVector(vector<vector<TYPE>> vecList)
{
    vector<TYPE> result;
    
    int n = 0;
    for (int i = 0; i < vecList.size(); i++)
    {
        n += vecList[i].size();
    }

    result.reserve(n);
    for (int i = 0; i < vecList.size(); i++)
    {
        result.insert(result.end(), vecList[i].begin(), vecList[i].end());
    }
    
    return result;
}

class SignalEvent
{
public:
    TLorentzVector higgs1, higgs2, hardJet;
    TLorentzVector b1, b2, other1, other2;
    vector<TLorentzVector> jetList;
    bool status;
    int nJet;
    double Ht, thrust, topness, dih_bbaa, deltaPhi, topness2;

    SignalEvent()
    {
        jetList = {};
    }
    ~SignalEvent()
    {

    }

    double diHiggsInvM(void);
};

struct Cut
{
    Cut()
    {
        deltaR_bb_max = 100.0;
        deltaR_bb_min = 0.0;
        deltaR_aa_max = 100.0;
        deltaR_aa_min = 0.0;
        deltaR_ab_max = 100.0;
        deltaR_ab_min = 0.0;
        delta_maa = 100000.0;
        delta_mbb = 100000.0;
        ptb = 0.0;
        ptj = 0.0;
        fakePhoton = fakebJet = 0;
        topness = 0;
        topness2 = 0;
        deltaPhi = HRatio = 0;
    }
    double deltaR_bb_max;
    double deltaR_bb_min;
    double deltaR_aa_max;
    double deltaR_aa_min;
    double deltaR_ab_max;
    double deltaR_ab_min;
    double delta_maa;
    double delta_mbb;
    double ptb;
    double ptj;
    double topness, topness2;
    int fakePhoton, fakebJet;
    double deltaPhi;
    double HRatio;
};




#endif