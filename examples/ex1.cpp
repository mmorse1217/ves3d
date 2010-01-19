#include "vesicle.h"
#include<iostream>

int main(int argc, char* argv[])
{
    vesicle<double> V(1,13);
    
    for(int ii=0;ii<V.p;++ii)
	std::cout<<V.posVec[ii]<<std::endl;

    return 0;
}
