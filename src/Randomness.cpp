#include <limits>
#include <random>
#include <cstring>

#include "Randomness.hpp" 
#include "Density.hpp"


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
  ++ctr.v[0]; 
  if(ctr.v[0] == std::numeric_limits<unsigned int>::max())
    {
      ++ctr.v[1]; 
      ctr.v[0] = 0;     
    }
}



int Randomness::drawIntegerClosed(int upperBound)
{
  assert(upperBound >= 0 ); 
  if(upperBound == 0)
    return 0 ; 
  else 
    {
      randCtr_t r = exa_rand(key, ctr); 
      ++ctr.v[0]; 
      return r.v[0] % (upperBound + 1 ) ; 
    }
}


int Randomness::drawIntegerOpen(int upperBound)
{
  // cout << "drawing from [0," << upperBound << ")" << endl; 
  assert(upperBound > 0); 
  return drawIntegerClosed(upperBound-1);   
}


nat Randomness::operator()() 
{
  auto seed = generateSeed(); 
  return seed.v[0]; 
}


double Randomness::drawRandDouble01()
{
  randCtr_t r = exa_rand(key, ctr); 
  ++ctr.v[0]; 

  return r123::u01<double>(r.v[1]);
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


double Randomness::drawMultiplier(double multiplier) // 
{

  double tmp =  exp(multiplier * (drawRandDouble01()  - 0.5)); 
  assert(tmp > 0.); 
  return tmp ;   
}


double Randomness::drawFromSlidingWindow(double value, double window)
{
  return value + window * (drawRandDouble01() - 0.5); 
}


void Randomness::drawPermutation( int* perm, int n)
{
  int i;
  int randomNumber;
  perm[0] = 0;
  
  for(i=1 ; i<n ; i++){
  
    randomNumber = drawIntegerOpen(i+1);
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
std::vector<double> Randomness::drawRandDirichlet( const std::vector<double> &alphas)
{
  // we have to set some lower bound. If the gamma distribution
  // proposes 0 (and it will), everythings falls apart if the
  // normalization rate (the last one) becomes too small.
  // 
  // The boundary is a rather small value anyway.
  static double lowerBound = 1e-32; 
  
  auto sum= double{0};
  auto result = std::vector<double>{}; 

  for(auto alpha : alphas)
    {
      double val = drawRandGamma(alpha, 1.0); 
      if(val <= lowerBound)
	val = lowerBound; 

      result.push_back(val); 
      sum += val;
    }

  for_each(begin(result), end(result), [&](double &v) {  v /= sum; }); 

  return result; 
}


double Randomness::drawRandGamma(double alpha, double beta)
{
  // a hack, until we have our own gamma 
  auto generator = std::default_random_engine{} ;
  generator.seed((*this)());   
  auto dist = std::gamma_distribution<double>(alpha,1./beta); // 1/beta necessary?
  return dist(generator); 
}


nat Randomness::drawGeometric(double prop)  
{
  auto gen = std::default_random_engine{} ; 
  gen.seed((*this)());
  auto dist = std::geometric_distribution<nat>(prop);
  return dist(gen); 
}


nat Randomness::drawBinomial(double prop, nat trials )
{
  auto dist = std::binomial_distribution<nat>(trials,prop);
  return dist(*this);
}


std::ostream& operator<<(std::ostream& out, const Randomness &rhs)
{
  out << "key={" << rhs.key.v[0] << "," << rhs.key.v[1] << "},ctr={" << rhs.ctr.v[0] << ","<< rhs.ctr.v[1] << "}"; 
  return out; 
} 



void Randomness::deserialize( std::istream &in )  
{ 
  key.v[0] = cRead<nat>(in); 
  key.v[1] = cRead<nat>(in) ; 
  ctr.v[0] = cRead<nat>(in); 
  ctr.v[1] = cRead<nat>(in); 
}
 
void Randomness::serialize( std::ostream &out)  const
{
  cWrite(out, key.v[0]); 
  cWrite(out, key.v[1]); 
  cWrite(out, ctr.v[0]); 
  cWrite(out, ctr.v[1]); 
} 


void Randomness::rebaseForGeneration(nat generation)
{
  ctr.v[1] = generation; 
  ctr.v[0] = 0; 
}



nat Randomness::getGeneration() const 
{
  return ctr.v[1]; 
}


void Randomness::setKey(randKey_t _key) 
{
  key = _key; 
}
