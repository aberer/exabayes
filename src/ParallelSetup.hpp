/** 
    @file ParallelSetup.hpp
    @brief handles all high-level parallel communication    
 */ 

#ifndef _PARALLEL_SETUP_H
#define _PARALLEL_SETUP_H


#include "config.h"

class SampleMaster;
class CoupledChains; 

#include <iostream>


#if HAVE_PLL == 0
#include "ParallelSetupExaml.hpp"
#else 
#include "ParallelSetupPll.hpp"
#endif
#endif
