#include "RandomVariable.hpp"

map<category_t, string> RandomVariable::nameMap = {
    { TOPOLOGY, "topo" } , 
    { BRANCH_LENGTHS, "bl" },
    { FREQUENCIES , "pi" } ,
    { SUBSTITUTION_RATES, "revMat"} , 
    { RATE_HETEROGENEITY , "shape" } , 
    { AA_MODEL, "aaModel" } 
  } ; 


ostream& operator<<(ostream &out, const RandomVariable& rhs)
{
  out << RandomVariable::nameMap[rhs.cat] ; 

  bool  isFirst = true; 
  out << "{" ; 
  for(auto v : rhs.partitions)
    {
      if (not isFirst) 
	out << "," ; 
      else 
	isFirst = false; 
      out << v; 
    }
  out << "} \twith prior " << rhs.getPrior();   
  return out;
}
