#include "config.h" 


#ifdef HAVE_AVX

#define  __AVX

#if HAVE_PLL == 1 
#include "pll/avxLikelihood-pll.c"
#else 
#include "examl/avxLikelihood.c"
#endif
#else 


int prototypeDummy()
{
  return 0; 
}

#endif
