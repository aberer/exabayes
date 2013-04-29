#include <math.h>
#include <inttypes.h>

#include <random>

#include "axml.h"
#include "main-common.h"
#include "bayes.h"
#include "globals.h"
#include "randomness.h"
#include "branch.h"
#include "adapters.h"

#include "TreeAln.hpp" 




#if 0 
static void inc_global()
{
  if(gAInfo.rGlobalCtr.v[0] == 2147483645) /* roughly what a int32_t can hold */
    {
      gAInfo.rGlobalCtr.v[1]++; 
      gAInfo.rGlobalCtr.v[0] = 0;       
    }
  else 
    gAInfo.rGlobalCtr.v[0]++; 
  
  assert(gAInfo.rGlobalCtr.v[1] < 2147483645); 
}


randCtr_t drawGlobalRandInt()
{
  randCtr_t result = exa_rand(gAInfo.rGlobalKey, gAInfo.rGlobalCtr); 
  inc_global();
  return result; 
}

int drawGlobalRandIntBound(int upperBound)
{
  randCtr_t r = drawGlobalRandInt();
  return r.v[0] % upperBound; 
}



double drawGlobalDouble01()
{
  randCtr_t result  = exa_rand(gAInfo.rGlobalKey, gAInfo.rGlobalCtr);   
  return u01_closed_open_32_53(result.v[1]) ; 
}

#endif

/* chain specific  */




/* uniform r in  [0,upperBound-1] */
// int drawRandInt(Chain *chain, int upperBound )
// {
//   chain->rCtr.v[0] =  chain->currentGeneration; 
//   randCtr_t r = exa_rand(chain->rKey, chain->rCtr); 
//   chain->rCtr.v[1]++;   
//   return r.v[0] % upperBound; 
// }


/* uniform r in [0,1) */







  
  
  







