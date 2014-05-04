#include "log_double.hpp"


log_double exponentiate(log_double val, double exponent) 
{
  val._val *= exponent; 
  return val; 
}


std::ostream& operator<<(std::ostream& out, const log_double& rhs )
{
  out << rhs._val ;
  return out; 
}
