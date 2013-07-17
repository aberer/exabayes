#ifndef _PARAMETER_CONTENT
#define _PARAMETER_CONTENT

#include <iostream>
#include <vector>

#include "Branch.hpp"
#include "Checkpointable.hpp"

class  ParameterContent : public Checkpointable
{
public:   
  std::vector<double> values; 
  std::vector<Branch> branches; 

  virtual void readFromCheckpoint( std::ifstream &in )  ; 
  virtual void writeToCheckpoint( std::ofstream &out) ;   

  // AA model? 

  friend std::ostream& operator<<(std::ostream& out,  const ParameterContent &rhs); 
}; 

#endif
