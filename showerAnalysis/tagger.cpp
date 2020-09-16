#include "iostream"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJetStructureBase.hh"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

using namespace fastjet;
using namespace std;
using namespace HepMC;


vector<PseudoJet> getFinalState(GenEvent *event)
{
    vector<PseudoJet> finalState;
    int n = 1;
    for (GenEvent::particle_iterator p = event->particles_begin() ; p != event->particles_end(); ++p)
    {
        if (!((*p)->end_vertex()) && ((*p)->status() == 1))
        {
            PseudoJet tmp = PseudoJet((*p)->momentum().x(), (*p)->momentum().y(), (*p)->momentum().z(), (*p)->momentum().e());
            tmp.set_user_index(n);
            finalState.push_back(tmp);
            n += 1;
        }
    }
    return finalState;
}

vector<PseudoJet> massDrop(vector<PseudoJet> inputList)
{
    vector<PseudoJet> outputList, fatjets;
    vector<PseudoJet>::iterator itOutputList;
    ClusterSequence *sequence = new ClusterSequence(inputList, JetDefinition(cambridge_algorithm, 0.4));
    fatjets = sorted_by_pt(sequence->inclusive_jets(40));
    //cout << "fatjet number" << " " << fatjets.size() << endl;
    //cout << "fatjet mass" << " " << fatjets[0].m();
    std::cout << inputList.size() << endl;
    int n = 0;
    for (int i = 0; i < fatjets.size(); i++)
    {
        n += fatjets[i].constituents().size();
    }
    std::cout << n << endl; 
    
    if (fatjets.size() > 0)
    {
        vector<PseudoJet> fatSubstructure;
        PseudoJet parent1, parent2;
        vector<PseudoJet> parentJet = {fatjets[0]};
        vector<PseudoJet>::iterator itParentJet;

        while (parentJet.size() > 0)
        {
            vector<PseudoJet> tmp;
            for (itParentJet = parentJet.begin(); itParentJet != parentJet.end(); ++itParentJet)
            {
                if (!(*itParentJet).has_parents(parent1, parent2)) continue;
                /* cout << parent1.m() << endl;
                cout << parent2.m() << endl; */
                if (parent1.m() > parent2.m() && parent1.m() > 0.8 * (*itParentJet).m())
                {
                    tmp.push_back(parent1);
                }
                else if (parent2.m() > parent1.m() && parent2.m() > 0.8 * (*itParentJet).m())
                {
                    tmp.push_back(parent2);
                }
                else
                {
                    tmp.push_back(parent1);
                    tmp.push_back(parent2);
                }
                for (int i = 0; i < tmp.size(); i++)
                {
                    if (tmp[i].m() < 30)
                    {
                        fatSubstructure.push_back(tmp[i]);
                        tmp.erase(tmp.begin() + i);
                    }
                }
            }
            parentJet = tmp;
            tmp.clear();
        }
        //cout << "fatjet substructure number" << " " << fatSubstructure.size() << endl;
        // remove fatjet constituents
        vector<PseudoJet> myInputList = inputList, fatConstituents = fatjets[0].constituents();
        vector<PseudoJet>::iterator itMyInput;
        vector<int> fatJetIndex;
        for (int i = 0; i < fatConstituents.size(); i++)
        {
            fatJetIndex.push_back(fatConstituents[i].user_index());
        }
        //cout << fatJetIndex.size() << endl; 
        //cout << myInputList.size() << endl;
        for (itMyInput = myInputList.begin(); itMyInput != myInputList.end(); ++itMyInput)
        {
            //cout << myInputList.size() << endl;
            for (int i = 0; i < fatJetIndex.size(); i++)
            {
                if ((*itMyInput).user_index() == fatJetIndex[i])
                {
                    myInputList.erase(itMyInput);
                    --itMyInput;
                    break;
                } 
            }
        }
        //cout << myInputList.size() << endl;
        fatJetIndex.clear();
        // re-cluster remmants
        ClusterSequence *remmantSequence;
        vector<PseudoJet> remmantList;
        double remmantPt = 30;
        remmantSequence = new ClusterSequence(myInputList, JetDefinition(antikt_algorithm, 0.4));
        remmantList = sorted_by_pt(remmantSequence->inclusive_jets(remmantPt));
        //cout << "rest jet number" << " " << remmantList.size() << endl;
        // write results to outputlist
        outputList = fatSubstructure;
        for (int i = 0; i < remmantList.size(); i++)
        {
            outputList.push_back(remmantList[i]);
        }
        for (itOutputList = outputList.begin(); itOutputList != outputList.end(); ++itOutputList)
        {
            if (itOutputList < outputList.begin() + fatSubstructure.size())
            {
                //cout << "True" << endl;
            }
            
        }
        
        
    }
    return outputList;
}


int main(int argc, char *argv[])
{
    IO_GenEvent ascii("/mnt/d/work/Hpair/events/hepmc/hhj1.hepmc", ios::in);
    GenEvent *event = ascii.read_next_event();
    //int n = 0;
    /* for (int i = 0; i < 3; i++)
    {
        ascii >> event;
    } */
    /* while (event)
    {
        massDrop(getFinalState(event));
        ascii >> event;
    } */
    massDrop(getFinalState(event));
    return 0;
}