/** 
    @file Sampling.hpp
    
    @brief represents various chains sampling the posterior probability space
    
    Despite of its modest name, this is in fact the master class.  
 */ 


#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <vector>

#include "CoupledChains.hpp"

using namespace std; 


class Sampling
{
public: 
  Sampling();
  ~Sampling();

private: 
  vector<CoupledChains> runs; 

  
};  


#endif
