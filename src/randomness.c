/**
   @file randomness.c
   
   @brief Wraps random number generation functions. 

   @notice Do not generate random numbers otherwise! 
 */ 


#include <math.h>
#include <assert.h>
#include <inttypes.h>


#include "common.h"
#include "axml.h"
#include "main-common.h"
#include "proposalStructs.h"
#include "globals.h"
#include "randomness.h"





/* 
   notice new convention: 
   
   * chain-specific random numbers (e.g., proposals, acceptance) are
     created from randCtr_t and randKey_t and => this way we easily
     can reset the rng

   * random numbers that are needed in a global context (e.g.,
     switching chains) are created from a global rng. (see global
     variables)
     
   * for the ctr for the chain-specific stuf: first int is the
     generation, second a ctr (since we may need more numbers for each
     step, starting from 0)

 */


/* here are a few wrapper functions: we will not be using the standard
   random number generator forever. */


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




void initLocalRng(state *theChain)
{
  /* init rng */
  theChain->rCtr.v[0] = 0; 
  theChain->rCtr.v[1] = 0; 
  randCtr_t r = drawGlobalRandInt();
  theChain->rKey.v[0] = r.v[0]; 
  theChain->rKey.v[1] = r.v[1]; 
  if(processID == 0)
    printf("initialized chain %d with seed %d,%d\n", theChain->id, theChain->rKey.v[0], theChain->rKey.v[1]); 
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



/* chain specific  */



/* uniform r in  [0,upperBound-1] */
int drawRandInt(state *chain, int upperBound )
{
  chain->rCtr.v[0] =  chain->currentGeneration; 
  randCtr_t r = exa_rand(chain->rKey, chain->rCtr); 
  chain->rCtr.v[1]++;   
  return r.v[0] % upperBound; 
}


/* uniform r in [0,1) */
double drawRandDouble01(state *chain)
{
  chain->rCtr.v[0] = chain->currentGeneration; 
  randCtr_t r = exa_rand(chain->rKey, chain->rCtr); 
  chain->rCtr.v[1]++; 
  return u01_closed_open_32_53(r.v[1]); 
}



double drawRandExp(state *chain, double lambda)
{  
  double r = drawRandDouble01(chain);   
  return -log(r )/ lambda; 
}


double drawRandBiUnif(state *chain, double x)
{
  double r = drawRandDouble01(chain) *  (2*x-x/2) + x / (3/2) ; 
  return r; 
}


//Given x this function randomly draws from [x/2,2*x] 
/* double drawRandBiUnif(double x) */
/* { */
/*   double r; */
/*   //r=(double)rand()/(double)RAND_MAX;   */
/*   //r=x/2+r*(3/2)*x;//=x/2+r*(2*x-x/2); */
  
/*   r=drawRandDouble(2*x-x/2)+x/2; */
/*   //r=drawRandDouble((3/2)*x-x/(3/2))+x/(3/2); */
  
/*   return r; */
/* } */


/**
   @brief draw r according to distribution given by weights. 

   NOTE sum of weights is not required to be 1.0
*/
int drawSampleProportionally(state *chain,  double *weights, int numWeight )
{
  double r = drawRandDouble01(chain);
 
  double sum=0.0;
  double lower_bound = 0.0;
  size_t i = 0;
  
  assert( numWeight > 0 );
  
  for( i = 0; i < numWeight ; ++i ) 
    {
      sum+=weights[i]; 
    }
  assert(sum>0);
  r=r*sum;
    
  for( i = 0; i < numWeight ; ++i ) 
    {
      double upper_bound = lower_bound + weights[i];
    
      if( r >= lower_bound && r < upper_bound ) 
	return i ;
    
      lower_bound = upper_bound; 
    }
  
  return i-1;
}

//get random permutation of [0,n-1]
void drawPermutation(state *chain, int* perm, int n)
{
  int i;
  int randomNumber;
  perm[0] = 0;
  
  for(i=1 ; i<n ; i++){
  
    randomNumber = drawRandInt(chain, i+1);
    // randomNumber=rand() % (i+1);
    // randomNumber=rand();

    if(randomNumber==i){
      perm[i]=i;
    }else{
      perm[i]=perm[randomNumber];
      perm[randomNumber]=i;
    }
  }
  
  /*for(i=0 ; i<n ; i++){
    printf("%d ",perm[i]);

    
    }
    printf("\n");
  */
} 


