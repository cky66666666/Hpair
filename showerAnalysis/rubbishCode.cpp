#include <vector>
#include "iostream"

using namespace std;

void test()
{
    vector<double> a = {1,2,3};
    vector<double>::iterator itA;
    for (itA = a.begin(); itA != a.end() && itA != a.end() + 1; ++itA)
    {
        cout << a.size() << endl;
        if (*itA == 3)
        {
            a.erase(itA);
        }
    }
    cout << a.size() << endl;
    //cout << *itA << endl;
}

int main()
{
    for (int i = 0; i < 3; i++)
    {
        test();
    }
       
}