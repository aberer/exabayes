#include "Parameters.hpp"

category_t RevMatParameter::cat = SUBSTITUTION_RATES; 
category_t FrequencyParameter::cat = FREQUENCIES; 
category_t RateHetParameter::cat = RATE_HETEROGENEITY; 

bool RevMatParameter::needsBLupdate = true; 
bool FrequencyParameter::needsBLupdate = true; 
bool RateHetParameter::needsBLupdate = false; 


// double RevMatParameter::initWeight  = 0.5 ; 
// double FrequencyParameter::initWeight  = 0.5 ; 
// double RateHetParameter::initWeight  = 1 ; 
