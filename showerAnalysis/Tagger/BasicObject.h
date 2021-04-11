#ifndef _BasicObject_
#define _BasicObject_

#include "vector"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

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
    int nJet;
    double Ht;

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
        deltaR_bb_max = 2.0;
        deltaR_bb_min = 0.4;
        deltaR_aa_max = 2.0;
        deltaR_aa_min = 0.4;
        deltaR_ab_max = 2.0;
        deltaR_ab_min = 0.4;
        delta_maa = 3.0;
        delta_mbb = 25.0;
        ptb = 40.0;
        ptj = 100.0;
        nljet = 1000.0;
        nnljet = 1000.0;
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
    double nljet;
    double nnljet;
};




#endif