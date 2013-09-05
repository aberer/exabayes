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
  assert(0); 

  // if(rhs.values.size() > 0)
  //   for_each(rhs.values.begin(), rhs.values.end(), [&](const double &d) {out << d << "," ; }) ; 
  // else if(rhs.branchLengths.size() > 0)
  //   for_each(rhs.branchLengths.begin(), rhs.branchLengths.end(), [&](const BranchLength &b) {out << b << "," ; }) ; 

  return out; 
} 
