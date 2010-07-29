#include <vector>
#include <sys/time.h>
#include <cmath>
#include "stokes.h"
#include <cassert>
#include <iostream>
#include <cstdlib>

#define get_seconds()   (gettimeofday(&tp, &tzp), \
    (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)

struct timeval  tp;
struct timezone tzp;

using namespace std;

int main(int argc, char ** argv)
{
  srand48(0*time(0));

  int m=(12+1)*2*12;
  int Nv = 100000;
  vector<float> src(3*Nv*m,1);
  vector<float> den(3*Nv*m,1);
  vector<float> pot(3*Nv*m,12345);
  vector<float> pot2(3*Nv*m,12345);
  vector<float> q(m,1);

  for (size_t i=0; i<src.size(); i++)
  {
    src[i] = i%5+1;
    den[i] = i%6+1;
  }


  double ss = get_seconds();
  DirectStokes(m, Nv, 0, 1, &q[0], &src[0], &src[0], &den[0], &pot[0]);
  cout<<"Stokes NO SSE:   "<<get_seconds()-ss<<" sec\n";


  ss = get_seconds();
  DirectStokesSSE(m, Nv, 0, 1, &q[0], &src[0], &src[0], &den[0], &pot2[0]);
  cout<<"Stokes WITH SSE: "<<get_seconds()-ss<<" sec\n";

  for (size_t i=0; i<pot.size(); i++)
  {
    float base = max( fabs (pot[i]),fabs (pot2[i]) );
    if (base>0)
      if ( fabs(pot[i]-pot2[i]) > 1e-5*base )
      {
	cout << "i= "<<i<<" nosse= "<<pot[i]<<" sse= "<<pot2[i]<<endl;
	abort();
      }
  }
  return 0;
}

