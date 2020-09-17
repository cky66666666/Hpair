#include "Tagger/BoostedHiggs.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "vector"

using namespace std;
using namespace fastjet;

class HiggsCandidate
{
public:
    HiggsCandidate(vector<PseudoJet> candJet, vector<PseudoJet> candConstituent)
    {
        this->candJet = candJet;
        this->candConstituent = candConstituent;
    }
    vector<PseudoJet> candJet, candConstituent;
};

BoostedHiggs::BoostedHiggs(vector<PseudoJet> finalState)
{
    this->finalState = finalState;
}

BoostedHiggs::~BoostedHiggs()
{

}

void BoostedHiggs::findFatJet()
{
    fatSequence = new ClusterSequence(this->finalState, JetDefinition(cambridge_algorithm, 1.5));
    fatJet= sorted_by_pt(fatSequence->inclusive_jets(150.0));

    vector<PseudoJet> mother, tmp;
    vector<PseudoJet>::iterator itMother;
    PseudoJet parent1, parent2;

    if (fatJet.size() > 0)
    {
        if (fatJet[0].m() > 110)
        {
            mother = {fatJet[0]};
        }
        else
        {
            fatSubStructure = {};
        }
    }
    else
    {
        fatSubStructure = {};
    }
    
    while (mother.size() > 0)
    {
        for (itMother = mother.begin(); itMother != mother.end(); ++itMother)
        {
            if (!(*itMother).has_parents(parent1, parent2)) continue;
            if (parent1.m() > parent2.m() && parent1.m() > 0.8 * (*itMother).m())
            {
                tmp.push_back(parent1);
            }
            else if (parent2.m() > parent1.m() && parent2.m() > 0.8 * (*itMother).m())
            {
                tmp.push_back(parent2);
            }
            else
            {
                tmp.push_back(parent1);
                tmp.push_back(parent2);
            }
        }
        mother.clear();
        for (int i = 0; i < tmp.size(); i++)
        {
            if (tmp[i].m() < 30)
            {
                fatSubStructure.push_back(tmp[i]);
            }
            else
            {
                mother.push_back(tmp[i]);
            }
            
        }
        tmp.clear();
    }
}