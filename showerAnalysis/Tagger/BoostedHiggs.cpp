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

BoostedHiggs::BoostedHiggs(vector<PseudoJet> finalState, vector<GenParticle*> parton, double fatJetMinPt = 150, double fatJetMinMass = 110, double fatJetR = 1.5, double fatSubStructureR = 0.3)
{
    this->finalState = finalState;
    this->fatJetMinPt = fatJetMinPt;
    this->fatJetR = fatJetR;
    this->fatJetMinMass = fatJetMinMass;
    this->fatSubStructureR = fatSubStructureR;
    this->parton = parton;
}

BoostedHiggs::~BoostedHiggs()
{

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

void BoostedHiggs::declusterFatJet()
{
    fatSequence = new ClusterSequence(finalState, JetDefinition(cambridge_algorithm, fatJetR));
    fatJet= sorted_by_pt(fatSequence->inclusive_jets(fatJetMinPt));

    vector<PseudoJet> mother, tmp;
    vector<PseudoJet>::iterator itMother;
    PseudoJet parent1, parent2;

    if (fatJet.size() > 0)
    {
        if (fatJet[0].m() > fatJetMinMass)
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

void BoostedHiggs::getHiggsCandidate()
{
    int nSubStructure = fatSubStructure.size();
    vector<PseudoJet> combiedJet, constituent1, constituent2, tmp;
    ClusterSequence *filter;
    double deltaR;
    higgsCandidateList = {};
    for (int i = 0; i < nSubStructure; i++)
    {
        for (int j = i + 1; j < nSubStructure; j++)
        {
            constituent1 = fatSubStructure[i].constituents();
            constituent2 = fatSubStructure[j].constituents();
            combiedJet.reserve(constituent1.size() + constituent2.size());
            combiedJet.insert(combiedJet.end(), constituent1.begin(), constituent1.end());
            combiedJet.insert(combiedJet.end(), constituent2.begin(), constituent2.end());

            deltaR = min(fatSubStructureR, deltaRCalc(fatSubStructure[i], fatSubStructure[j]) / 2);
            filter = new ClusterSequence(combiedJet, JetDefinition(cambridge_algorithm, deltaR));
            tmp = sorted_by_pt(filter->inclusive_jets());
            
            if (tmp.size() >= 3)
            {
                higgsCandidateList.push_back(HiggsCandidate({tmp[0], tmp[1], tmp[2]}, {tmp[0].constituents(), tmp[1].constituents(), tmp[2].constituents()}));
            }
        }
    }
}

void BoostedHiggs::findBoostedHiggs()
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
                deltaInvMass = abs((candidate[0] + candidate[1]).m() - 125);
                boostedHiggs = {candidate[0], candidate[1]};
                higgsConstituent = {};
                higgsConstituent.reserve(constituent[0].size() + constituent[1].size());
                higgsConstituent.insert(higgsConstituent.end(), constituent[0].begin(), constituent[0].end());
                higgsConstituent.insert(higgsConstituent.end(), constituent[1].begin(), constituent[1].end());
            }
            
        }
        
    }
    
}

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

void BoostedHiggs::process()
{
    declusterFatJet();
    getHiggsCandidate();
    findBoostedHiggs();
    getRemmant();
}