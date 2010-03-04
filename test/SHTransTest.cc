#include<iostream>
#include "DeviceCPU.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SHTrans.h"

int main(int argc, char ** argv)
{
    int p=12;
    int n_surfs = 1;
    
    DeviceCPU<float> cpu1, cpu2;
    Scalars<float> sc(cpu1, p, n_surfs);

    SHTrans<float> diff(cpu1, p, n_surfs);
    
    diff.AllDerivatives(sc,sc,sc,sc,sc,sc);
    diff.FirstDerivatives(sc,sc,sc);

    return 0;
}
