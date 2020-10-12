
#include "iostream"
#include "vector"

using namespace std;

int main()
{
    vector<double> test;
    for (int i = 0; i < 100000; i++)
    {
        test.push_back(i);
    }
    for (int i = 0; i < test.size(); i++)
    {
        cout << test[i];
    }
    
    
    return 0;
}