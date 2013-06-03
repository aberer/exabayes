#include "Parameters.hpp"


category_t RevMatParameter::cat = SUBSTITUTION_RATES; 
category_t FrequencyParameter::cat = FREQUENCIES; 
category_t RateHetParameter::cat = RATE_HETEROGENEITY; 


bool RevMatParameter::needsFcUpdate = true; 
bool FrequencyParameter::needsFcUpdate = true; 
bool RateHetParameter::needsFcUpdate = false; // correct? 
