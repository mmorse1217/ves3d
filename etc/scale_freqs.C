#include "scale_freqs.h"
#include <cassert>
#include <algorithm>

void scale_freqs(int p, int num_vesicles, float * inputs, float * alphas, float * outputs )
{
  float * inp_deb = inputs;
  float * out_deb = outputs;
  float * alphas_deb = alphas;


  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int leg_order = p+1;
  for (int v=0; v<num_vesicles; v++)
    for (int i=0; i<leg_order; i++)
      *(outputs++) = *(inputs++) * alphas[i];
  alphas += leg_order;
  leg_order--;

  // process remaining frequencies except the last cosine
  for (; leg_order>1; leg_order--) 
  {
    // first process cosine
    for (int v=0; v<num_vesicles; v++)
      for (int i=0; i<leg_order; i++)
	*(outputs++) = *(inputs++) *alphas[i];
    alphas += leg_order;

    // then process sine
    for (int v=0; v<num_vesicles; v++)
      for (int i=0; i<leg_order; i++)
	*(outputs++) = *(inputs++) *alphas[i];
    alphas += leg_order;
  }

  // process last cosine
  for (int v=0; v<num_vesicles; v++)
    *(outputs++) = *(inputs++) * alphas[0];
  alphas += leg_order;
  leg_order--;

  assert (leg_order == 0);
  assert(inputs-inp_deb == num_vesicles*p*(p+2));
  assert(outputs-out_deb == num_vesicles*p*(p+2));
  assert(alphas-alphas_deb == p*(p+2));
}

void resample(int p, int num_vesicles, int q, float * inputs, float * outputs)
{
  float * inp_deb = inputs;
  float * out_deb = outputs;



  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int leg_order = p+1;
  int new_leg_order = q+1;
  int min_leg_order = std::min(leg_order, new_leg_order);

  for (int v=0; v<num_vesicles; v++)
  {
    for(int i=0; i<min_leg_order; i++)
      *(outputs++) = *(inputs++);
    for (int i=leg_order; i<new_leg_order; i++)
      *(outputs++) = 0;
    if (leg_order > new_leg_order)
      inputs += leg_order - new_leg_order;
  }
  leg_order--;
  new_leg_order--;
  min_leg_order--;

  // process remaining frequencies except the last cosine
  for (; min_leg_order>1; min_leg_order--,leg_order--,new_leg_order--) 
  {
    // first process cosine
    for (int v=0; v<num_vesicles; v++)
    {
      for(int i=0; i<min_leg_order; i++)
	*(outputs++) = *(inputs++);
      for (int i=leg_order; i<new_leg_order; i++)
	*(outputs++) = 0;
      if (leg_order > new_leg_order)
	inputs += leg_order - new_leg_order;
    }

    // then process sine
    for (int v=0; v<num_vesicles; v++)
    {
      for(int i=0; i<min_leg_order; i++)
	*(outputs++) = *(inputs++);
      for (int i=leg_order; i<new_leg_order; i++)
	*(outputs++) = 0;
      if (leg_order > new_leg_order)
	inputs += leg_order - new_leg_order;
    }
  }

  // process last cosine
  for (int v=0; v<num_vesicles; v++)
  {
    for(int i=0; i<min_leg_order; i++)
      *(outputs++) = *(inputs++);
    for (int i=leg_order; i<new_leg_order; i++)
      *(outputs++) = 0;
    if (leg_order > new_leg_order)
      inputs += leg_order - new_leg_order;
  }

  leg_order--;
  new_leg_order--;
  min_leg_order--;

  // assert (leg_order == 0);
 
  // if q>p all remaining coefs should be zero
  float * output_end = out_deb+num_vesicles*q*(q+2);
  assert(outputs<=output_end);

  while (outputs<output_end)
    *(outputs++) = 0;

  if (p<=q)
    assert(inputs-inp_deb == num_vesicles*p*(p+2));
  else
    assert(inputs-inp_deb < num_vesicles*p*(p+2));
    
  assert(outputs-out_deb == num_vesicles*q*(q+2));
}

