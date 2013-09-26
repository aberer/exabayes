#include "AbstractParameter.hpp" 
#include "Category.hpp"


AbstractParameter::AbstractParameter(Category cat, nat id, nat _idOfMyKind)
  : id(id)
  , idOfMyKind(_idOfMyKind)
  , cat(cat) 
  , printToParamFile(true)
{
  
}


std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs)
{
  return rhs->printShort(out);
}


std::ostream&  AbstractParameter::printShort(std::ostream& out) const
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

