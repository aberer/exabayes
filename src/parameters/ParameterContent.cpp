#include "ParameterContent.hpp"
#include <limits>
#include <iomanip>
#include <algorithm> 


void ParameterContent::readFromCheckpoint( std::istream &in )  
{
  nat found = 0; 
  if(values.size() > 0 )
    ++found; 
  if(branchLengths.size() > 0)
    ++found; 
  if(topology.size() > 0)
    ++found; 
  assert(found == 1 ); 

  for(auto &v : values)
    v = cRead<double>(in); 

  for(auto &b : branchLengths)
    b.readFromCheckpoint(in);

  for(auto &b : topology)
    b.readFromCheckpoint(in); 
} 


void ParameterContent::writeToCheckpoint( std::ostream &out) const 
{
  nat found = 0; 
  if(values.size() > 0 )
    ++found; 
  if(branchLengths.size() > 0)
    ++found; 
  if(topology.size() > 0)
    ++found; 
  assert(found == 1 ); 

  for(auto &v : values)
    cWrite(out, v); 

  for(auto &b : branchLengths)    
    b.writeToCheckpoint(out);

  for(auto &b : topology)
    b.writeToCheckpoint(out); 
}   



std::ostream& operator<<(std::ostream& out, const ParameterContent &rhs)
{
  if(rhs.values.size() > 0)
    out << rhs.values; 
  else if(rhs.branchLengths.size( )> 0)
    out << rhs.branchLengths; 
  else if(rhs.topology.size() > 0)
    out << rhs.topology ; 
  else 
    assert(0); 

  return out; 
} 
