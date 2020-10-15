#ifndef _BAChannel_
#define _BAChannel_

#include "vector"
#include "TLorentzVector.h"

namespace fastjet
{
class PseudoJet;
class ClusterSequence;
}

class GenParticle;
class Photon;
class HiggsCand;

using namespace std;
using namespace fastjet;

class BAChannel
{
private:
    vector<PseudoJet> finalState;
    vector<GenParticle*> parton;
    vector<TLorentzVector> photon, electron, muon, jet;
    vector<HiggsCand> higgsCandListA, higgsCandListB;
    
    HiggsCand higgsSelector(vector<HiggsCand>);
    int flavourAssociation(TLorentzVector);
    bool preselect();

public:
    TLorentzVector higgsFromA, higgsFromB, hardJet;
    bool status;

    BAChannel();
    ~BAChannel();

    void init(vector<PseudoJet>, vector<GenParticle*>, vector<TLorentzVector>, vector<TLorentzVector>, vector<TLorentzVector>);
    void preprocess();
    void find2AHiggs();
    void find2BHiggs();
    void find2BHiggsHard();
    void process();
    void finish();

};

#endif