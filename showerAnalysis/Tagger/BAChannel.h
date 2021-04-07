#ifndef _BAChannel_
#define _BAChannel_

#include "vector"
#include "TLorentzVector.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "classes/DelphesClasses.h"
#include "vector"
#include "BoostedHiggs.h"
#include "BasicObject.h"

using namespace std;
using namespace fastjet;

class HiggsCand
{
public:
    HiggsCand()
    {

    }
    ~HiggsCand()
    {

    }
    TLorentzVector cand1;
    TLorentzVector cand2;
    int index1, index2;
};

class BAChannel
{
private:
    vector<PseudoJet> finalState;
    vector<GenParticle*> parton;
    vector<TLorentzVector> photon, electron, muon, bJet, lightJet;
    vector<Jet*> delphesJet;
    vector<HiggsCand> higgsCandListA, higgsCandListB;
    vector<TLorentzVector> photonPair, bPair;
    Cut cut;
    
    HiggsCand higgsSelector(vector<HiggsCand>);
    int flavourAssociation(TLorentzVector);
    vector<TLorentzVector> ptSort(vector<TLorentzVector>);
    bool trigger();

public:
    SignalEvent signal;
    bool status;

    BAChannel();
    ~BAChannel();

    void init(vector<PseudoJet>, vector<GenParticle*>, vector<TLorentzVector>, vector<TLorentzVector>, vector<TLorentzVector>, vector<Jet*>, Cut cut);
    void preprocess();
    void selPhotonPair();
    void selBPair();
    void find2BHiggsHard();
    void process();
    void finish();

};

#endif