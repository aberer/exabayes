#include "Randomness.hpp" 

#include "TreeAln.hpp"
#include "branch.h"

#include <limits>

#include "densities.h"

// TODO proper AND carefull make-over of randomness


Randomness::Randomness(int seed)
{
  key.v[0] = seed; 
  key.v[1] = 0; 
  ctr.v[0] = 0; 
  ctr.v[1] = 0; 
}

randCtr_t Randomness::generateSeed()
{
  randCtr_t r =  exa_rand(ctr, key);   
  incrementNoLimit();
  return r; 
}




void Randomness::incrementNoLimit()
{
  ctr.v[0]++; 
  if(ctr.v[0] == std::numeric_limits<unsigned int>::max())
    {
      ctr.v[1]++; 
      ctr.v[0] = 0;     
    }
}


int Randomness::drawRandInt(int upperBound )
{
  randCtr_t r = exa_rand(key, ctr); 
  ctr.v[1]++;   
  return r.v[0] % upperBound; 
}


double Randomness::drawRandDouble01()
{
  randCtr_t r = exa_rand(key, ctr); 
  ctr.v[1]++; 
  return u01_closed_open_32_53(r.v[1]); 
}


double Randomness::drawRandExp(double lambda)
{  
  double r = drawRandDouble01();   
  return -log(r )/ lambda; 
}



double Randomness::drawRandBiUnif(double x)
{
  double r = drawRandDouble01() *  (2*x-x/2) + x / (3/2) ; 
  return r; 
}




/** @brief gets a multiplier for updating a parameter or branch length */
double Randomness::drawMultiplier(double multiplier) // 
{

  double tmp =  exp(multiplier * (drawRandDouble01()  - 0.5)); 
  assert(tmp > 0.); 
  return tmp ;   
}



double Randomness::drawFromSlidingWindow(double param, double window)
{
  double upper = param + (window / 2 ) ,
    lower = param - (window / 2 ); 
  
  double r = drawRandDouble01(); 
  return lower + r * (upper - lower)  ; /* TODO correct?  */
}




/**
   @brief draws a subtree uniformly. 

   We have to account for tips again. 
   
   The thisNode contains the root of the subtree.
*/ 
branch Randomness::drawSubtreeUniform(TreeAln &traln)
{
  tree *tr = traln.getTr();   
  while(TRUE)
    {
      branch b = drawBranchUniform(traln); 
      if(isTipBranch(b,tr->mxtips))
	{
	  if(isTip(b.thisNode, tr->mxtips))
	    b = invertBranch(b); 
	  if(drawRandDouble01() < 0.5)
	    return b; 
	}
      else 
	{
	  return drawRandDouble01() < 0.5 ? b  : invertBranch(b); 	  
	}
    }
}



/**
   @brief draws a branch with uniform probability.
   
   We have to treat inner and outer branches separatedly.
 */
branch Randomness::drawBranchUniform(TreeAln &traln)
{
  tree *tr = traln.getTr(); 

  boolean accept = FALSE; 
  int randId = 0; 
  while(NOT accept)
    {
      randId = drawRandInt(2 * tr->mxtips - 2 ) + 1;
      assert(randId > 0); 
      double r = drawRandDouble01(); 
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
      int r = drawRandInt(2); 
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
   @brief draws an inner branch 
*/ 
branch Randomness::drawInnerBranchUniform(TreeAln &traln)
{
  tree *tr = traln.getTr(); 
  boolean accepted = FALSE; 
  branch b; 
  do 
    {
      b = drawBranchUniform(traln);
      accepted = NOT isTip(b.thisNode, tr->mxtips ) && NOT isTip(b.thatNode, tr->mxtips);       
    }while(NOT accepted); 

  return b; 
}



//get random permutation of [0,n-1]
void Randomness::drawPermutation( int* perm, int n)
{
  int i;
  int randomNumber;
  perm[0] = 0;
  
  for(i=1 ; i<n ; i++){
  
    randomNumber = drawRandInt(i+1);
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
   @brief draw r according to distribution given by weights. 

   NOTE sum of weights is not required to be 1.0
*/
int Randomness::drawSampleProportionally( double *weights, int numWeight )
{
  double r = drawRandDouble01();
 
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











//This function should be called if the alphas for the dirichlet distribution are given
void Randomness::drawRandDirichlet( double* results, double* alphas,  int length)
{
  double sum=0;
  for(int i=0; i< length;i++)
    {
      results[i]=drawRandGamma(alphas[i], 1.0);
      sum+=results[i];
    }
  for(int i=0; i< length;i++)
    {
      results[i]=results[i]/sum;
    }
}



//This function should be called if the expected values for the dirichlet distribution are given
void Randomness::drawDirichletExpected(double* results, double* mean,double scale, int length)
{
  double alphas[length];
  double originalSum=0;

  for(int i=0; i< length;i++)
    {
      originalSum+=mean[i];
      alphas[i]=mean[i]*scale;
    }
  
  drawRandDirichlet( results, alphas, length);

  for(int i=0; i< length;i++)
    {
      results[i]=results[i]*originalSum;    
    }
  
}







// Gamma(alpha, beta) sampling
double Randomness::drawRandGamma(double alpha, double beta)
{

  double alpha_min=0.0001;
  if(alpha<alpha_min)
    alpha=alpha_min;
  
  if(beta<0)
    beta=0;
  
  double gamma=0;
   
  int escape=0, limit=2000;
  bool alert=false;
   
  if (alpha < 1) {
    double r;
    boolean done=false;
    double b=1+1/alpha;
    while(!done) {
      r=drawRandDouble01();
      if (r>1/b) {
	gamma=-log((b-r)/alpha);
	if (drawRandDouble01()<=pow(gamma,alpha-1)) 
	  done=true;
      }
      else {
	gamma=pow(r,1/alpha);
	if (drawRandDouble01()<=exp(-gamma)) 
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
    gamma = -log (drawRandDouble01());
  } else {    
    double y = -log (drawRandDouble01());
    while (drawRandDouble01() > pow (y * exp (1 - y), alpha - 1)){
      y = -log (drawRandDouble01());
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



ostream& operator<<(ostream& out, const Randomness &rhs)
{
  out << "key={" << rhs.key.v[0] << "," << rhs.key.v[1] << "},ctr={" << rhs.ctr.v[0] << ","<< rhs.ctr.v[1] << "}"; 
  return out; 
} 

