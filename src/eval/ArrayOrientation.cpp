#include "ArrayOrientation.hpp"
#include <cassert>
#include "Branch.hpp"


ArrayOrientation::ArrayOrientation(const TreeAln &traln)
  : orientation(traln.getNumberOfPartitions())
{
  for(nat i = 0 ; i < traln.getNumberOfPartitions(); ++i)
    orientation[i] = std::vector<nat>(traln.getNumberOfInnerNodes(),  INVALID); 
}



void ArrayOrientation::setOrientationForAllPartitions(nat id ,nat value)  
{
  for(nat i= 0; i < orientation.size() ; ++i)
    orientation[i][id] = value; 
}



void ArrayOrientation::setPartitionInvalid(nat part)  
{
  for(auto &elem  : orientation[part])
    elem = INVALID; 
} 


std::ostream& operator<<(std::ostream& out, const ArrayOrientation &rhs)
{
  nat ctr = 0; 
  for(auto &elem : rhs.orientation)
    {
      out << "[" << ctr << "] " ; 
      
      nat ctr2 = rhs.orientation[0].size() + 3 ;
      for(auto &subelem : elem ) 
	{
	  if(subelem == 0)
	    out << ctr2 <<  "->INVALID, " ; 
	  else 
	    out << ctr2 << "->" << subelem << ", " ;
	  ++ctr2; 
	}
      ++ctr ; 
      out << std::endl; 
    }
  return out; 
}


void ArrayOrientation::setInvalid(nat part, nat id)
{
  orientation[part][id] = INVALID; 
} 


nat ArrayOrientation::getOrientation(nat part, nat id ) const 
{
  auto value = orientation.at(part).at(id); 
  return value; 
}

