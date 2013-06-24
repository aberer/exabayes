#include "RandomVariable.hpp"



ostream& operator<<(ostream &out, const RandomVariable& rhs)
{
  out << CategoryFuns::getShortName(rhs.cat); 

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


ostream&  RandomVariable::printShort(ostream& out)
{
  out << CategoryFuns::getShortName(cat) << "{" ; 
  bool isFirst= true; 
  for(auto v : partitions)
    {
      if(not isFirst)
	out << ","; 
      else 
	isFirst = false; 
      out  << v ; 
    }
  out << "}";     
  return out; 
}


RandomVariablePtr RandomVariable::clone() const 
{
  RandomVariablePtr r(new RandomVariable(cat, id)) ; 
  r->setPrior(prior) ; 
  for(auto v : partitions)
    r->addPartition(v); 
  return r;  
}
