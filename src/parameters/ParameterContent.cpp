#include "ParameterContent.hpp"
#include <limits>
#include <iomanip>
#include <algorithm> 


void ParameterContent::readFromCheckpoint( std::ifstream &in )  
{
  assert(not (values.size() > 0 && branches.size() > 0) ); 
  
  for(auto &v : values)
    {
      in >> v ; 
      readDelimiter(in);
    }

  for(auto &b : branches)
    b.readFromCheckpoint(in);
} 


void ParameterContent::writeToCheckpoint( std::ofstream &out) 
{
  assert(not (values.size() > 0 && branches.size() > 0) ); 
  out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 2 );
  
  for(auto &v : values)
    out << v << DELIM; 

  for(auto &b :branches)
    b.writeToCheckpoint(out);
}   



std::ostream& operator<<(std::ostream& out, const ParameterContent &rhs)
{
  if(rhs.values.size() > 0)
    for_each(rhs.values.begin(), rhs.values.end(), [&](const double &d) {out << d << "," ; }) ; 
  else if(rhs.branches.size() > 0)
    for_each(rhs.branches.begin(), rhs.branches.end(), [&](const Branch &b) {out << b << "," ; }) ; 

  return out; 
} 
