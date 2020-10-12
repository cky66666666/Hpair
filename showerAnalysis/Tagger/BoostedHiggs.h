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
    BoostedHiggs();
    ~BoostedHiggs();

    void init(std::vector<fastjet::PseudoJet>, std::vector<GenParticle*>, double, double, double, double);
    void clear();

    void declusterFatJet();
    void getHiggsCandidate();
    std::vector<fastjet::PseudoJet> findBoostedHiggs(int);
    void getRemmant();
    
    void process(int);

    double deltaRCalc(fastjet::PseudoJet, fastjet::PseudoJet);
    int flavourAssociation(fastjet::PseudoJet);

    std::vector<fastjet::PseudoJet> remmant, boostedHiggs;
    double fatJetMinPt, fatJetMinMass, fatJetR, fatSubStructureR;

private:
    fastjet::ClusterSequence *fatSequence;
    std::vector<fastjet::PseudoJet> finalState, fatSubStructure, higgsConstituent, fatJet;
    std::vector<GenParticle*> parton;
    std::vector<HiggsCandidate> higgsCandidateList;
    
    bool isHiggs(std::vector<fastjet::PseudoJet>);
    std::vector<fastjet::PseudoJet> combineJet(std::vector<std::vector<fastjet::PseudoJet>>);
};


#endif