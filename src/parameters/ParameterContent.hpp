#ifndef _PARAMETER_CONTENT
#define _PARAMETER_CONTENT

#include <iostream>
#include <vector>

#include "Branch.hpp"
#include "Serializable.hpp"

 

class  ParameterContent : public Serializable
{
public:   
  std::vector<double> values; 
  std::vector<BranchPlain> topology;
  std::vector<BranchLength> branchLengths; 

  virtual void deserialize( std::istream &in )  ; 
  virtual void serialize( std::ostream &out)  const;   

  // AA model? 

  friend std::ostream& operator<<(std::ostream& out,  const ParameterContent &rhs); 
}; 

#endif
