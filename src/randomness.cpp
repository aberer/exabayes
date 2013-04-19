


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
#if 0 
  if(processID == 0)
    printf("initialized chain %d with seed %d,%d\n", theChain->id, theChain->rKey.v[0], theChain->rKey.v[1]); 
#endif
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


void normalize(double * vector, int length, double normalizingConstant)
{
  double sum=0;
  for(int i=0; i<length; i++)
    sum+=vector[i];
  
  for(int i=0; i<length; i++)
    vector[i]=vector[i]*normalizingConstant/sum;
  
}

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


double gammaFunction(double alpha)
{
  
double  lgam = lgamma(alpha);
 return signgam*exp(lgam); 
}

double betaFunction(double *alpha, int length)
{
 double beta=1.0;
 double sum=0;
  for(int i=0; i<length;i++)
  {
   beta*=gammaFunction(alpha[i]); 
  sum+=alpha[i];    
  }
  beta=beta/(gammaFunction(sum));
 
 return beta;
}

//density for dirichlet distribution with parameters "alphas" at point "values".
double densityDirichlet(double *values, double *alphas, int length)
{
  double density=1;
  double normValues[length];
  
  for(int i=0; i<length; i++)
  normValues[i]=values[i];
  
  density=density/betaFunction(alphas, length);
  
  normalize(normValues, length, 1);
  
  for(int i=0; i<length; i++)
    density=density*pow(normValues[i],alphas[i]-1);
  
 return density; 
}

  // Gamma(alpha, beta) sampling
double drawRandGamma(state *chain, double alpha, double beta)
{

 double alpha_min=0.0001;
  if(alpha<alpha_min)
    alpha=alpha_min;
  
  if(beta<0)
    beta=0;
  
   double gamma=0;
   
   int escape=0, limit=2000;
   boolean alert=true;
   
  if (alpha < 1) {
      double r;
      boolean done=false;
      double b=1+1/alpha;
      while(!done) {
        r=drawRandDouble01(chain);
     if (r>1/b) {
          gamma=-log((b-r)/alpha);
          if (drawRandDouble01(chain)<=pow(gamma,alpha-1)) 
	    done=true;
        }
	else {
          gamma=pow(r,1/alpha);
          if (drawRandDouble01(chain)<=exp(-gamma)) 
	    done=true;
        }
        
        
        
	    escape++;
	    if(alert && escape>limit)
	    {
	    printf("r>1/b==%.2f>%.2f\n",r,1/b);
	    printf("FREQ_MIN==%f alpha==%f gamma==%.2f\n",FREQ_MIN, alpha, gamma);
	    alert=false; 
	    //assert(0);
	    }
	
	
	
       }
    } 
    else if (alpha == 1) {//GAMMA(1,1)=EXP(1)
      gamma = -log (drawRandDouble01(chain));
    } else {    
      double y = -log (drawRandDouble01(chain));
      while (drawRandDouble01(chain) > pow (y * exp (1 - y), alpha - 1)){
        y = -log (drawRandDouble01(chain));
	escape++;
	if(alert && escape>limit){
	  printf("in second\n");
	 alert=false; 
	}
      }
      gamma = alpha * y;
    }
    return beta*gamma;//scale from GAMMA(alpha,1) to GAMMA(alpha,beta)
  }  
  
  
  //This function should be called if the alphas for the dirichlet distribution are given
void drawRandDirichlet(state *chain, double* results, double* alphas,  int length)
{
  double sum=0;
  for(int i=0; i< length;i++)
  {
    results[i]=drawRandGamma(chain, alphas[i], 1.0);
    sum+=results[i];
  }
   for(int i=0; i< length;i++)
  {
    results[i]=results[i]/sum;
  }
}

  //This function should be called if the expected values for the dirichlet distribution are given
void drawDirichletExpected(state *chain, double* results, double* mean,double scale, int length)
{
  double alphas[length];
  double originalSum=0;

     for(int i=0; i< length;i++)
  {
    originalSum+=mean[i];
    alphas[i]=mean[i]*scale;
  }
  
  drawRandDirichlet(chain, results, alphas, length);

  for(int i=0; i< length;i++)
  {
  results[i]=results[i]*originalSum;    
  }
  
}


/**
   @brief draw r according to distribution given by weights. 

   NOTE sum of weights is not required to be 1.0
*/
int drawSampleProportionally(state *chain,  double *weights, int numWeight )
{
  double r = drawRandDouble01(chain);
 
  double sum=0.0;
  double lower_bound = 0.0;
  int i = 0; 
  
  assert( numWeight > 0 );
  
  for(  i = 0; i < numWeight ; ++i ) 
    {
      sum+=weights[i]; 
    }
  assert(sum>0);
  r=r*sum;
    
  for( int i = 0; i < numWeight ; ++i ) 
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





/**
    @brief draws an inner branch 
 */ 
branch drawInnerBranchUniform(state *chain)
{
  tree *tr = chain->traln->getTr(); 
  boolean accepted = FALSE; 
  branch b; 
  do 
    {
      b = drawBranchUniform(chain);
      accepted = NOT isTip(b.thisNode, tr->mxtips ) && NOT isTip(b.thatNode, tr->mxtips);       
    }while(NOT accepted); 

  return b; 
}




/**
   @brief draws a subtree uniformly. 

   We have to account for tips again. 
   
   The thisNode contains the root of the subtree.
*/ 
branch drawSubtreeUniform(state *chain)
{
  tree *tr = chain->traln->getTr();   
  while(TRUE)
    {
      branch b = drawBranchUniform(chain); 
      if(isTipBranch(b,tr->mxtips))
	{
	  if(isTip(b.thisNode, tr->mxtips))
	    b = invertBranch(b); 
	  if(drawRandDouble01(chain) < 0.5)
	    return b; 
	}
      else 
	{
	  return drawRandDouble01(chain) < 0.5 ? b  : invertBranch(b); 	  
	}
    }
}


/**
   @brief draws a branch with uniform probability.
   
   We have to treat inner and outer branches separatedly.
 */
branch drawBranchUniform(state *chain)
{
  tree *tr = chain->traln->getTr(); 

  boolean accept = FALSE; 
  int randId = 0; 
  while(NOT accept)
    {
      randId = drawRandInt(chain, 2 * tr->mxtips - 2 ) + 1;
      assert(randId > 0); 
      double r = drawRandDouble01(chain); 
      if(isTip(randId, tr->mxtips) )
	accept = r < 0.25 ; 
      else 
	accept = r <= 0.75; 	
    }

  branch result; 
  result.thisNode = randId; 
  nodeptr p = tr->nodep[randId]; 
  if(isTip(randId, tr->mxtips))
    result.thatNode = p->back->number; 
  else 
    {
      int r = drawRandInt(chain,2); 
      switch(r)
	{
	case 0 : 
	  result.thatNode = p->back->number; 
	  break; 
	case 1 : 
	  result.thatNode = p->next->back->number; 
	  break; 
	case 2: 
	  result.thatNode = p->next->next->back->number; 
	  break; 
	default: assert(0); 
	}
    }
  
  return result; 
}


/**
   @brief gets a multiplier for updating a parameter or branch length 
 */
double drawMultiplier(state *chain, double multiplier)
{
  double tmp =  exp(multiplier * (drawRandDouble01(chain)  - 0.5)); 
  assert(tmp > 0.); 
  return tmp ;   
}


double drawFromSlidingWindow(state *chain, double param, double window)
{
  double upper = param + (window / 2 ) ,
    lower = param - (window / 2 ); 
  
  double r = drawRandDouble01(chain); 
  return lower + r * (upper - lower)  ; /* TODO correct?  */
}
