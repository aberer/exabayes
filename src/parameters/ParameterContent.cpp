#include "ParameterContent.hpp"
#include <limits>
#include <iomanip>


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


void ParameterContent::writeToCheckpoint( std::ofstream &out) const
{
  assert(not (values.size() > 0 && branches.size() > 0) ); 
  out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 2 );
  
  for(auto &v : values)
    out << v << DELIM; 

  for(auto &b :branches)
    b.writeToCheckpoint(out);
}   
