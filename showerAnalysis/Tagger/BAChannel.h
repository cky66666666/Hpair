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

using namespace std;
using namespace fastjet;

struct HiggsCand;
class BAChannel
{
private:
    vector<PseudoJet> finalState;
    vector<GenParticle*> parton;
    vector<TLorentzVector> photon, electron, muon, jet, bJet;
    vector<HiggsCand> higgsCandListA, higgsCandListB;
    
    int flavourAssociation(TLorentzVector);

public:
    TLorentzVector higgsFromA, higgsFromB, hardJet;
    bool status;

    BAChannel();
    ~BAChannel();

    void init(vector<PseudoJet>, vector<GenParticle*>, vector<TLorentzVector>, vector<TLorentzVector>, vector<TLorentzVector>);
    void preprocess();
    void find2AHiggs();
    void find2BHiggs();
    void process();

};

BAChannel::BAChannel()
{
}

BAChannel::~BAChannel()
{
}


#endif