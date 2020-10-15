#include "Tagger/BoostedHiggs.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TLorentzVector.h"

#include "vector"
#include "iostream"

#include "classes/DelphesClasses.h"

using namespace std;
using namespace fastjet;

class HiggsCandidate
{
public:
    HiggsCandidate(vector<PseudoJet> candJet, vector<vector<PseudoJet>> candConstituent)
    {
        this->candJet = candJet;
        this->candConstituent = candConstituent;
    }
    vector<vector<PseudoJet>> candConstituent;
    vector<PseudoJet> candJet;
};

BoostedHiggs::BoostedHiggs()
{
    
}

BoostedHiggs::~BoostedHiggs()
{
    
}

void BoostedHiggs::init(vector<PseudoJet> finalState, vector<GenParticle*> parton, vector<PseudoJet> lepHiggs, double fatJetMinPt = 150, double fatJetMinMass = 110, double fatJetR = 1.5, double fatSubStructureR = 0.3)
{
    this->finalState = finalState;
    this->fatJetMinPt = fatJetMinPt;
    this->fatJetR = fatJetR;
    this->fatJetMinMass = fatJetMinMass;
    this->fatSubStructureR = fatSubStructureR;
    this->parton = parton;
    this->lepHiggs = lepHiggs;

    boostedHiggs = {};
    higgsConstituent = {};
}

void BoostedHiggs::clear()
{
    higgsConstituent.clear();
    fatJet.clear();
    boostedHiggs.clear();
    remmant.clear();
    finalState.clear();
    fatSubStructure.clear();
    parton.clear();
    higgsCandidateList.clear();
    delete fatSequence;
}

double BoostedHiggs::deltaRCalc(PseudoJet jet1, PseudoJet jet2)
{
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(jet1.px(), jet1.py(), jet1.pz(), jet1.e());
    p2.SetPxPyPzE(jet2.px(), jet2.py(), jet2.pz(), jet2.e());
    return p1.DeltaR(p2);
}

int BoostedHiggs::flavourAssociation(PseudoJet jet)
{
    int flavour = -1;
    double deltaR = 0.5;
    TLorentzVector pJet;
    pJet.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.e());
    for (int i = 0; i < parton.size(); i++)
    {
        int pid = abs(parton[i]->PID);
        if (pid == 21) pid = 0;
        if (pJet.DeltaR(parton[i]->P4()) < deltaR && pid > flavour)
        {
            flavour = pid;
        }
    }
    return flavour;
}

bool BoostedHiggs::isHiggs(vector<PseudoJet> higgsCandidate)
{
    // higgsCandidate should be sorted according to pt
    double mH = 125;
    vector<PseudoJet> bJet, candidate;
    if (higgsCandidate.size() < 3) return false;
    candidate = {higgsCandidate[0], higgsCandidate[1], higgsCandidate[2]};

    for (int i = 0; i < candidate.size(); i++)
    {
        if (flavourAssociation(candidate[i]) == 5)
        {
            bJet.push_back(candidate[i]);
        }
    }

    if (bJet.size() != 2) return false;
    if (/* bJet[0].delta_R(bJet[1]) >= 0.4 && */ abs((candidate[0] + candidate[1] + candidate[2]).m() - mH) < 20)
    {
        return true;
    }
    else
    {
        return false;
    }
    

    /* if (abs(flavourAssociation(higgsCandidate[0])) == 5 && abs(flavourAssociation(higgsCandidate[1])) == 5 && abs((higgsCandidate[0] + higgsCandidate[1] + higgsCandidate[2]).m() - mH) < 20)
    {
        return true;
    }
    else
    {
        return false;
    } */
    
}

vector<PseudoJet> BoostedHiggs::combineJet(vector<vector<PseudoJet>> jetList)
{
    int n = jetList.size(), jetSize = 0;
    vector<PseudoJet> combinedJet;
    for (int i = 0; i < n; i++)
    {
        jetSize += jetList[i].size();
    }
    combinedJet.reserve(jetSize);
    for (int i = 0; i < n; i++)
    {
        combinedJet.insert(combinedJet.end(), jetList[i].begin(), jetList[i].end());
    }
    return combinedJet;
}

void BoostedHiggs::declusterFatJet()
{
    fatSequence = new ClusterSequence(finalState, JetDefinition(cambridge_algorithm, fatJetR));
    fatJet = sorted_by_pt(fatSequence->inclusive_jets(fatJetMinPt));

    vector<PseudoJet> mother = {}, tmp;
    vector<PseudoJet>::iterator itMother;
    PseudoJet parent1, parent2;

    if (lepHiggs.size() == 2)
    {
        for (int i = 0; i < fatJet.size(); i++)
        {
            if (fatJet[i].m() > fatJetMinMass && fatJet[i].delta_R(lepHiggs[0]) > 0.4 && fatJet[i].delta_R(lepHiggs[1]) > 0.4)
            {
                mother = {fatJet[i]};
                break;
            }
        }
    }
    else
    {
        for (int i = 0; i < fatJet.size(); i++)
        {
            if (fatJet[i].m() > fatJetMinMass)
            {
                mother = {fatJet[i]};
                break;
            }
        }
    }
    
    //cout << mother.size() << endl;

    fatSubStructure = {};
    
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

vector<PseudoJet> BoostedHiggs::findBoostedHiggs(int nHiggs)
{
    int nSubStructure = fatSubStructure.size();
    vector<PseudoJet> combinedJet, tmp, result;
    double deltaR;

    result = {};
    
    for (int i = 0; i < fatSubStructure.size(); i++)
    {
        for (int j = i + 1; j < fatSubStructure.size(); j++)
        {
            combinedJet = combineJet({fatSubStructure[i].constituents(), fatSubStructure[j].constituents()});
            deltaR = min(fatSubStructureR, deltaRCalc(fatSubStructure[i], fatSubStructure[j]) / 2);
            ClusterSequence filter = ClusterSequence(combinedJet, JetDefinition(cambridge_algorithm, deltaR));
            tmp = sorted_by_pt(filter.inclusive_jets());
            
            if (isHiggs(tmp))
            {
                result.push_back(tmp[0] + tmp[1] + tmp[2]);
                higgsConstituent = combineJet({higgsConstituent, tmp[0].constituents(), tmp[1].constituents(), tmp[2].constituents()});
                fatSubStructure.erase(fatSubStructure.begin() + j);
                fatSubStructure.erase(fatSubStructure.begin() + i);
                i--;
                break;
            }
        }
        if (result.size() == nHiggs)
        {
            break;
        }
        
    }
    return result;
}

/* void BoostedHiggs::findBoostedHiggs()
{
    double deltaInvMass = 1000, mH = 125;
    vector<PseudoJet> candidate;
    vector<vector<PseudoJet>> constituent;
    boostedHiggs = {};
    for (int i = 0; i < higgsCandidateList.size(); i++)
    {
        double tmp = 1000;
        for (int j = 0; j < higgsCandidateList[i].candJet.size(); j++)
        {
            for (int k = 0; k < higgsCandidateList[i].candJet.size(); k++)
            {
                if (abs((higgsCandidateList[i].candJet[j] + higgsCandidateList[i].candJet[k]).m() - mH) < tmp)
                {
                    tmp = abs((higgsCandidateList[i].candJet[j] + higgsCandidateList[i].candJet[k]).m() - mH);
                    candidate = {higgsCandidateList[i].candJet[j], higgsCandidateList[i].candJet[k]};
                    constituent = {higgsCandidateList[i].candConstituent[j], higgsCandidateList[i].candConstituent[k]};
                }
                
            }
            
        }
        if (candidate.size() < 2) continue;
        if (flavourAssociation(candidate[0]) == 5 && flavourAssociation(candidate[1]) == 5)
        {
            if (abs((candidate[0] + candidate[1]).m() - mH) < deltaInvMass)
            {
                deltaInvMass = abs((candidate[0] + candidate[1]).m() - mH);
                boostedHiggs = {candidate[0], candidate[1]};
                higgsConstituent = {};
                higgsConstituent.reserve(constituent[0].size() + constituent[1].size());
                higgsConstituent.insert(higgsConstituent.end(), constituent[0].begin(), constituent[0].end());
                higgsConstituent.insert(higgsConstituent.end(), constituent[1].begin(), constituent[1].end());
            }
            
        }
        
    }
    
} */

void BoostedHiggs::getRemmant()
{
    remmant = finalState;
    int n = higgsConstituent.size();
    vector<PseudoJet>::iterator itRemmant;
    for (itRemmant = remmant.begin(); itRemmant != remmant.end(); itRemmant++)
    {
        for (int i = 0; i < n; i++)
        {
            if ((*itRemmant).user_index() == higgsConstituent[i].user_index())
            {
                remmant.erase(itRemmant);
                itRemmant--;
                break;
            }
        }
    }
}

void BoostedHiggs::process(int type)
{
    //type=1: single boosted higgs
    //type=2: double boosted higgs
    //type=3: two collinear higgs
    if (type == 1)
    {
        declusterFatJet();
        boostedHiggs = findBoostedHiggs(1);
        getRemmant();
    }
    else if (type == 2)
    {
        vector<PseudoJet> tmp1, tmp2;
        
        declusterFatJet();
        tmp1 = findBoostedHiggs(1);
        getRemmant();
        higgsConstituent = {};

        finalState = remmant;
        declusterFatJet();
        tmp2 = findBoostedHiggs(1);
        getRemmant();

        boostedHiggs = combineJet({tmp1, tmp2});
    }
    else if (type == 3)
    {
        declusterFatJet();
        boostedHiggs = findBoostedHiggs(2);
        getRemmant();
    }
    else
    {
        cout << "invalid input" << endl;
    }
}