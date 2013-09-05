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
  std::vector<BranchPlain> topology;
  std::vector<BranchLength> branchLengths; 

  virtual void readFromCheckpoint( std::istream &in )  ; 
  virtual void writeToCheckpoint( std::ostream &out)  const;   

  // AA model? 

  friend std::ostream& operator<<(std::ostream& out,  const ParameterContent &rhs); 
}; 

#endif
