#ifndef _SWAP_MATRIX_H 
#define _SWAP_MATRIX_H 

#include <vector>
#include <iostream>

#include "common.h"
#include "SuccessCounter.hpp"

class SwapMatrix
{
public: 
  SwapMatrix(nat numChains); 
  void update(nat a, nat b, bool acc); 
  const SuccessCounter& getCounter(nat a, nat b ) const; 

private: 
  nat mapToIndex(nat a, nat b) const ;

private:  
  std::vector<SuccessCounter> matrix; 
  nat numEntries; 


  friend  std::ostream& operator<<(std::ostream &out, SwapMatrix& rhs ) ; 
}; 


#endif
