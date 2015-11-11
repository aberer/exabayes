#include "time.hpp"


CLOCK::system_clock::time_point getTimePoint()
{
  return CLOCK::system_clock::now() ;
} 

double getDuration(CLOCK::system_clock::time_point tp )
{
  return CLOCK::duration_cast<CLOCK::duration<double> > (CLOCK::system_clock::now()- tp   ).count(); 
} 
