#include<iostream>
#include<exception>
#include "SHScalars.h"
#include "SHVectors.h"
#include "SphHarm.h"

int main(int argc, char ** argv)
{
  int p=12;
  int num_ves = 1;
  
  SHScalars<float> vecIn(p,num_ves);
  SphHarm<float> diff(p,num_ves);
  
  diff.AllDerivatives(vecIn,vecIn,vecIn,vecIn,vecIn,vecIn);
  diff.FirstDerivatives(vecIn,vecIn,vecIn);
}
