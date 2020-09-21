#ifndef _BOOSTEDHIGGS_
#define _BOOSTEDHIGGS_

#include "vector"

namespace fastjet
{
class PseudoJet;
class ClusterSequence;
}

class HiggsCandidate;
class TLorentzVector;
class GenParticle;


class BoostedHiggs
{
public:
    BoostedHiggs(std::vector<fastjet::PseudoJet>, std::vector<GenParticle*>, double, double, double, double);
    ~BoostedHiggs();

    void declusterFatJet();
    void getHiggsCandidate();
    void findBoostedHiggs();
    void getRemmant();
    void process();

    double deltaRCalc(fastjet::PseudoJet, fastjet::PseudoJet);
    int flavourAssociation(fastjet::PseudoJet);

    std::vector<fastjet::PseudoJet> higgsConstituent, fatJet, remmant, boostedHiggs;
    double fatJetMinPt, fatJetMinMass, fatJetR, fatSubStructureR;

private:
    fastjet::ClusterSequence *fatSequence;
    std::vector<fastjet::PseudoJet> finalState, fatSubStructure;
    std::vector<GenParticle*> parton;
    std::vector<HiggsCandidate> higgsCandidateList;
};


#endif