#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>




/* here are a few wrapper functions: we will not be using the standard
   random number generator forever. */


void initRNG(uint64_t seedHere)
{
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

