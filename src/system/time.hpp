#ifndef _TIME_EXA_HPP
#define _TIME_EXA_HPP

#include <chrono> 
#include <ratio>
#define CLOCK std::chrono  

CLOCK::system_clock::time_point getTimePoint(); 
double getDuration(CLOCK::system_clock::time_point tp ); 
#endif
