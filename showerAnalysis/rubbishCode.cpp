#include <vector>
#include "iostream"
#include "TH1D.h"
#include "fastjet/ClusterSequence.hh"

using namespace std;


int main()
{
    vector<double> test = {};
    TH1D *hist = new TH1D("a", "a", 50, 0, 10000);
    fastjet::ClusterSequence *sequence;
    int n = 0;

    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10000; j++)
        {
            hist->Fill(j);
        }
        
    }

    cout << n << endl;
}