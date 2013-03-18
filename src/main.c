/**
   @file main.c

   @brief A delegator file that contains the main function for builds
   with either ExaML or the PLL
 */ 

#include "config.h"

#ifdef HAVE_AVX
#define __AVX
#endif

#if HAVE_PLL == 1 
#include "src/main-pll_cond.c"
#else 
#include "src/main-examl_cond.c" 
#endif
