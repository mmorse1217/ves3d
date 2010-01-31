#include<iostream>
#include<exception>
#include "SHScalars.h"
#include "SHVectors.h"
#include "SphHarm.h"

using namespace std;

int main(int argc, char* argv[])
{
    SHVectors<double> vecIn(4,1);
    SphHarm<double> mamad;

    mamad.Derivatives(&vecIn,&vecIn,&vecIn,&vecIn,&vecIn,&vecIn);
}
