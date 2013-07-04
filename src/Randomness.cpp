#include <limits>

#include <cstring>

#include "Randomness.hpp" 
#include "densities.h"

// TODO proper AND carefull make-over of randomness


Randomness::Randomness(randCtr_t seed)
{
  key.v[0] = seed.v[0]; 
  key.v[1] = seed.v[1]; 
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




/** 
    @brief draw integer uniformly from [0,n]
 */ 
int Randomness::drawIntegerClosed(int upperBound)
{
  assert(upperBound >= 0 ); 
  if(upperBound == 0)
    return 0 ; 
  else 
    {
      randCtr_t r = exa_rand(key, ctr); 
      ctr.v[1]++;
      return r.v[0] % (upperBound + 1 ) ; 
    }
}


/** 
    @brief for transition: draw integer uniformly from [0,n)
 */ 
int Randomness::drawIntegerOpen(int upperBound)
{
  // cout << "drawing from [0," << upperBound << ")" << endl; 
  assert(upperBound > 0); 
  return drawIntegerClosed(upperBound-1);   
}



/** 
    @brief draw integer uniformly from [0,upperBound)
 */ 
int Randomness::drawRandInt(int upperBound )
{
  randCtr_t r = exa_rand(key, ctr); 
  ctr.v[1]++;   
  return r.v[0] % upperBound; 
}


nat Randomness::operator()() 
{
  randCtr_t r = exa_rand(key,ctr); 
  ctr.v[1]++;
  return r.v[0];  
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


// alternative slider as used in mrBayes
double Randomness::drawFromSlidingWindow(double value, double window)
{
  return value + window * (drawRandDouble01() - 0.5); 
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
void Randomness::drawRandDirichlet( std::vector<double> &results, const std::vector<double> &alphas)
{
  double sum=0;
  results.resize(alphas.size()); 
  for(nat i=0; i< alphas.size();i++)
    {
      results[i] = drawRandGamma(alphas[i], 1.0) ;
      sum += results[i];
    }
  for(nat i=0; i< alphas.size();i++)
    results[i] /= sum;
}



//This function should be called if the expected values for the dirichlet distribution are given
void Randomness::drawDirichletExpected(std::vector<double> &results, const std::vector<double> &mean,double scale)
{
  std::vector<double> alphas; 
  double originalSum=0;

  for(nat i=0; i< mean.size();i++)
    {
      originalSum+=mean[i];
      alphas.push_back( mean[i]*scale ) ;
    }
  
  drawRandDirichlet( results, alphas);

  for(nat i=0; i< alphas.size() ;i++)
    results[i]=results[i]*originalSum;      
}




// TODO : this is a bracen copy from mrBayes code 
#include <random>
#define ETA (1E-30)
double Randomness::drawRandGamma(double alpha, double beta)
{
#if 1
  // a hack, until we have our own gamma 
  std::default_random_engine generator;
  generator.seed((*this)());   
  std::gamma_distribution<double> dist(alpha,1/beta);

  return dist(generator); 
#else 


  double r= 0.0;
	
  if (alpha <= 0.0)    
    puts ("Gamma parameter less than zero\n");
  else if (alpha < 1.0)  
    {
      double r, x=0.0, small=1e-37, w;
      static double a, p, uf, ss=10.0, d;
	
      if (fabs(alpha-ss)>ETA) /* s != ss */ 
	{
	  a  = 1.0 - alpha;
	  p  = a / (a + alpha * exp(-a));
	  uf = p * pow(small / a, alpha);
	  d  = a * log(a);
	  ss = alpha;
	}
      for (;;) 
	{
	  cout << "stuck in first" << endl; 
	  r = drawRandDouble01();	  
	  if (r > p)        
	    x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;
	  else if (r>uf)  
	    x = a * pow(r / p, 1.0 / alpha), w = x;
	  else            
	    {
	      x  = 0; 
	      break; 
	    }
	  r = drawRandDouble01();
	  if (1.0 - r <= w && r > 0.0)
	    if (r*(w + 1.0) >= 1.0 || -log(r) <= w)  
	      continue;
	  break;
	}

      r = x ; 
    }
  else if (alpha > 1.0)  
    {
      double r , d, f, g, x;
      static double	b, h, ss=0.0;

      if (fabs(alpha-ss)>ETA) /* s != ss */
	{
	  b  = alpha - 1.0;
	  h  = sqrt(3.0 * alpha - 0.75);
	  ss = alpha;
	}
      for (;;) 
	{
	  r = drawRandDouble01();
	  g = r - r * r;
	  f = (r - 0.5) * h / sqrt(g);
	  x = b + f;
	  if (x <= 0.0) 
	    continue;
	  r = drawRandDouble01();
	  d = 64 * r * r * g * g * g;
	  if ( d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))  
	    break;
	}
		
      r = x ; 
    }
  else    
    r -= log(drawRandDouble01());
		
  return (r / beta);
#endif
}



// // is there something wrong with the code below? quarantined it for
// // now...
 
// // Gamma(alpha, beta) sampling
// double Randomness::drawRandGamma(double alpha, double beta)
// {
//   double alpha_min=0.0001;
//   if(alpha<alpha_min)
//     alpha=alpha_min;
  
//   if(beta<0)
//     beta=0;
  
//   double gamma=0;
   
//   int escape=0, limit=2000;
//   bool alert=false;
   
//   if (alpha < 1) {
//     double r;
//     boolean done=false;
//     double b=1+1/alpha;
//     while(!done) {
//       r=drawRandDouble01();
//       if (r>1/b) {
// 	gamma=-log((b-r)/alpha);
// 	if (drawRandDouble01()<=pow(gamma,alpha-1)) 
// 	  done=true;
//       }
//       else {
// 	gamma=pow(r,1/alpha);
// 	if (drawRandDouble01()<=exp(-gamma)) 
// 	  done=true;
//       }
        
        
        
//       escape++;
//       if(alert && escape>limit)
// 	{
// 	  printf("r>1/b==%.2f>%.2f\n",r,1/b);
// 	  printf("FREQ_MIN==%f alpha==%f gamma==%.2f\n",FREQ_MIN, alpha, gamma);
// 	  alert=false; 
// 	  //assert(0);
// 	}	
//     }
//   } 
//   else if (alpha == 1) {//GAMMA(1,1)=EXP(1)
//     gamma = -log (drawRandDouble01());
//   } else {    
//     double y = -log (drawRandDouble01());
//     while (drawRandDouble01() > pow (y * exp (1 - y), alpha - 1)){
//       y = -log (drawRandDouble01());
//       escape++;
//       if(alert && escape>limit){
// 	printf("in second\n");
// 	alert=false; 
//       }
//     }
//     gamma = alpha * y;
//   }
//   return beta*gamma;//scale from GAMMA(alpha,1) to GAMMA(alpha,beta)
// }



std::ostream& operator<<(std::ostream& out, const Randomness &rhs)
{
  out << "key={" << rhs.key.v[0] << "," << rhs.key.v[1] << "},ctr={" << rhs.ctr.v[0] << ","<< rhs.ctr.v[1] << "}"; 
  return out; 
} 
