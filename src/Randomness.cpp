#include "Randomness.hpp" 

#include <limits>

// TODO proper AND carefull make-over of randomness


Randomness::Randomness(int seed)
{
  key.v[0] = seed; 
  key.v[1] = 0; 
}

randCtr_t Randomness::generateSeed()
{
  randCtr_t r =  exa_rand(ctr, key);   
  incrementNoLimit();
  return r; 
}




void Randomness::incrementNoLimit()
{
  ctr.v[0]++; 
  if(ctr.v[0] == std::numeric_limits<unsigned int>::max())
    {
      ctr.v[1]++; 
      ctr.v[0] = 0;     
    }
}
