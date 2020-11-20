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
    int nJet;
    double Ht;

    SignalEvent()
    {

    }
    ~SignalEvent()
    {

    }

    double diHiggsInvM(void);
};




#endif