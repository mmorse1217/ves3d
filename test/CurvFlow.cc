#include<iostream>
#include "Vesicle.h"

using namespace std;

int main(int argc, char* argv[])
{

    Vesicle<double> S(12,1);
    *S.kappa_ = 1;
    //Set X for S
    SHVectors<double> bend_force(12,1);
    int m = 100;
    double ts = .1;

    for(int tt=0;tt<m;++tt)
    {
        AxPy(&S.k_, &S.normal_, 0.0, &bend_force);
        //Add pressure to preserve volume
        AxPy(ts,&bend_force, &S.x_, &S.x_);
    }

    return 0;
}
