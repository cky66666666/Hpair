#include <vector>
#include "iostream"
#include <omp.h>

using namespace std;


int main()
{
    vector<double> test = {};

    omp_set_num_threads(8);

    int n = 0;

    for (int i = 0; i < 10; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < 10000; j++)
        {
        //#pragma omp critical
            n += 1;
        }
        
    }

    cout << n << endl;
}