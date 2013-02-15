#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>

#include "common.h"
#include "axml.h"
#include "main-common.h"




/* here are a few wrapper functions: we will not be using the standard
   random number generator forever. */


void initRNG(uint64_t seedHere)
{
  PRINT("initializing RNG with seed %d\n", seedHere);   
  srand(seedHere);
}





/* uniform r in  [0,upperBound-1] */
int drawRandInt(int upperBound)
{
  return rand() % upperBound; 
}



/* uniform r in [0,1] */
double drawRandDouble()
{
  return (double)rand()/(double)RAND_MAX; 
}

double drawRandExp(double l)
{
  double r;
  r=(double)rand()/(double)RAND_MAX;
  r=-log(r)/l;
  return r; 
}

//Given x this function randomly draws from [x/2,2*x] 
double drawRandBiUnif(double x)
{
  double r;
  r=(double)rand()/(double)RAND_MAX;
  /*
   double s;
  s=(double)rand()/(double)RAND_MAX;
  
  if(r<1/3)//choose left area //if(r<1/3) yields exactly uniform on [x/2,2*x]
  {
   s=x/2+s*(x/2);//=x/2+s*(x-x/2) \in [x/2,x]
  }else{//choose right area
   s=x+s*x;//=x+s*(2*x-x) \in [x,2x]    
  }
  */
  r=x/2+r*(3/2)*x;//=x/2+r*(2*x-x/2);
  return r;
}


/* draw r according to distribution given by weights. NOTE sum of weights is not required to be 1.0*/
int drawSampleProportionally( double *weights, int numWeight )
{
  double r = drawRandDouble();
  
  double sum=0.0;
  float lower_bound = 0.0;
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
      float upper_bound = lower_bound + weights[i];
    
      if( r >= lower_bound && r < upper_bound ) 
	return i ;
    
      lower_bound = upper_bound; 
    }
  
  return i-1;
}

