#ifndef _PARAMETER_CONTENT
#define _PARAMETER_CONTENT

#include <vector>

#include "Branch.hpp"
#include "Checkpointable.hpp"

class  ParameterContent : public Checkpointable
{
public:   
  std::vector<double> values; 
  std::vector<Branch> branches; 

  virtual void readFromCheckpoint( std::ifstream &in )  ; 
  virtual void writeToCheckpoint( std::ofstream &out) const;   

  // AA model? 
}; 

#endif
