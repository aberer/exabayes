#ifndef _SWAP_MATRIX_H 
#define _SWAP_MATRIX_H 

#include <vector>
#include <iostream>

#include "common.h"
#include "SuccessCounter.hpp"

#include "Serializable.hpp"

class SwapMatrix : public Serializable
{
public: 
  SwapMatrix(nat numChains); 
  SwapMatrix(const SwapMatrix& rhs); 
  SwapMatrix operator=( SwapMatrix rhs ) ; 

  void update(nat a, nat b, bool acc); 
  const SuccessCounter& getCounter(nat a, nat b ) const; 

  virtual void deserialize( std::istream &in )   ; 
  virtual void serialize( std::ostream &out) const;   

  std::vector<SuccessCounter> getMatrix() const {return matrix; }

  SwapMatrix operator+(const SwapMatrix& rhs) const; 

  friend void swap(SwapMatrix& a, SwapMatrix &b); 

private: 
  nat mapToIndex(nat a, nat b) const ;

private:  
  std::vector<SuccessCounter> matrix; 
  nat numEntries; 


  friend  std::ostream& operator<<(std::ostream &out, SwapMatrix& rhs ) ; 
}; 


#endif
