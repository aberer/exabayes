#include "RandomVariable.hpp"



#if 0 
std::ostream& operator<<(std::ostream &out, const RandomVariable& rhs)
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


std::ostream&  RandomVariable::printShort(std::ostream& out)
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


#endif
