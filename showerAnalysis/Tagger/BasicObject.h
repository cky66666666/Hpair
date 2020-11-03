#ifndef _BasicObject_
#define _BasicObject_

#include "vector"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

using namespace std;

class SignalEvent
{
public:
    TLorentzVector higgs1, higgs2, hardJet;

    double diHiggsInvM(void);


};




#endif