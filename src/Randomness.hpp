#ifndef __RANDOMNESS_H
#define __RANDOMNESS_H


#include "rng.h"

class Randomness
{
public: 
  Randomness(int seed);
  randCtr_t generateSeed();   
  void rebase(int num){ ctr.v[1] = num; ctr.v[0] = 0; }
  

private:   
  randCtr_t ctr; 
  randKey_t key; 

  void incrementNoLimit();



}; 

#endif
