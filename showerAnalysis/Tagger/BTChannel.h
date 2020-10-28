#ifndef _BTChannel_
#define _BTChannel_

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "classes/DelphesClasses.h"
#include "TLorentzVector.h"

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TRandom2.h"

#include "vector"
#include "iostream"

using namespace fastjet;
using namespace std;

class MyEvent
{
public:
    MyEvent()
    {

    }
    ~MyEvent()
    {

    }
    TLorentzVector tau, antitau, b, antib, pMis, hardJet;
};

class BTChannel
{
private:
    vector<TLorentzVector> lepTau;
    vector<Jet*> tauJet, jet, bJet;
    vector<Electron*> electron;
    vector<Muon*> muon;
    TLorentzVector met;

    bool trigger(), selector();
    void findInvMom();
    double nutriMomFrac(TLorentzVector);
    double mT2();
    void consEvent();

    
public:
    BTChannel();
    ~BTChannel();

    bool status;
    double ditauInvM, dihiggsInvM;
    MyEvent event;

    void init(vector<Jet*>, vector<Electron*>, vector<Muon*>, TLorentzVector);
    double mTMax(const double*);
    void process();
    void finish();
};

#endif