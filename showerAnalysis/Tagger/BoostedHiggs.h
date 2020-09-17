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


class BoostedHiggs
{
public:
    BoostedHiggs(std::vector<fastjet::PseudoJet>);
    ~BoostedHiggs();

    void findFatJet();

    std::vector<fastjet::PseudoJet> constituent, boostedHiggs, fatJet;
    TLorentzVector *momentum;

private:


    fastjet::ClusterSequence *fatSequence;
    std::vector<fastjet::PseudoJet> finalState, fatSubStructure;
};


#endif