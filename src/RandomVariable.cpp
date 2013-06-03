#include "RandomVariable.hpp"

ostream& operator<<(ostream &out, const RandomVariable& rhs)
{
  map<category_t, string> nameMap = {
    { TOPOLOGY, "topo" } , 
    { BRANCH_LENGTHS, "bl" },
    { FREQUENCIES , "pi" } ,
    { SUBSTITUTION_RATES, "revMat"} , 
    { RATE_HETEROGENEITY , "shape" } , 
    { AA_MODEL, "aaModel" } 
  } ; 

  out << nameMap[rhs.cat] ; 

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
