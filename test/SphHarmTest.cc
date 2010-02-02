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
  SphHarm<float> mamad(p,num_ves);
  
  mamad.Derivatives(&vecIn,&vecIn,&vecIn,&vecIn,&vecIn,&vecIn);
}
