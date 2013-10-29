#include "ParameterContent.hpp"
#include <limits>
#include <iomanip>
#include <algorithm> 





void ParameterContent::deserialize( std::istream &in )  
{
  for(auto &v : values)
    v = cRead<double>(in); 

  for(auto &b : branchLengths)
    b.deserialize(in);

  for(auto &b : topology)
    b.deserialize(in); 
  
  for(auto &v :  protModel)
    v = ProtModel(cRead<int>(in));
} 


void ParameterContent::serialize( std::ostream &out) const 
{
  for(auto &v : values)
    cWrite(out, v); 

  for(auto &b : branchLengths)    
    b.serialize(out);

  for(auto &b : topology)
    b.serialize(out); 

  for(auto &v : protModel)
    {
      auto tmp = int(v);
      cWrite<int>(out,tmp); 
    }
}   



std::ostream& operator<<(std::ostream& out, const ParameterContent &rhs)
{
  if(rhs.values.size() > 0)
    out << rhs.values; 
  else if(rhs.branchLengths.size( )> 0)
    out << rhs.branchLengths; 
  else if(rhs.topology.size() > 0)
    out << rhs.topology ; 
  else if(rhs.protModel.size()  > 0)
    out << rhs.protModel; 
  else     
    assert(0); 

  return out; 
} 
