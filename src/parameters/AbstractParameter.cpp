#include "AbstractParameter.hpp" 
#include "Category.hpp"


std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs)
{
  out << CategoryFuns::getShortName(rhs->cat); 

  bool  isFirst = true; 
  out << "{" ; 
  for(auto v : rhs->partitions)
    {
      if (not isFirst) 
	out << "," ; 
      else 
	isFirst = false; 
      out << v; 
    }
  auto p = rhs->getPrior(); 
  out << "} \twith prior " << p ; 
  return out;
}


std::ostream&  AbstractParameter::printShort(std::ostream& out)
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
