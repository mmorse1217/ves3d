#include<iostream>
#include<exception>
#include "Vesicle.h"

using namespace std;

int main(int argc, char* argv[])
{

    Vesicle<double> S(4,5);

    *S.kappa_ = 1;
    cout<<*S.kappa_<<endl;
    return 0;
}
