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


log_double operator/(  const log_double& lhs , const log_double &rhs )
{
  return lhs * exponentiate(rhs,-1); 
}

log_double operator*( const log_double& lhs, const log_double &rhs)
{
#if 0 
  auto result = lhs; 
  result._error = lhs._error + rhs._error; 

  auto y = rhs._val - result._error; 
  auto t = result._val + y ; 
  result._error = (t - result._val)  -  y ; 
  result._val = t; 
  return result;  
#else 
  auto result = lhs ; 
  result._val += rhs._val; 
  return result; 
#endif 

}



log_double log_double::fromAbs(double val)
{
  if(  val <= 0 )
    {
      tout << "error: tried to create a log-double of " <<  val << std::endl; 
      assert(val > 0 ); 
    }
    
  auto result = log_double();
  result._val = log(val); 
  return result; 
}


log_double log_double::fromLog(double val )
{
  auto result = log_double();
  result._val = val; 
  return result; 
}

log_double log_double::lowest()
{
  auto result = log_double::fromAbs(1.); 
  result._val = std::numeric_limits<double>::lowest();
  return result; 
}
