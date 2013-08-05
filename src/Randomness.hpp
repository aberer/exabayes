#ifndef ___RANDOMNESS_H
#define ___RANDOMNESS_H

#include <iostream>
#include <vector>

#include "axml.h"

#include "GlobalVariables.hpp"
#include "Checkpointable.hpp" 

class TreeAln; 

/* TODO is this correct? */
#ifndef UINT64_C
#define UINT64_C
#endif

/* TODO
   
   we could mess around here a lot with 64-bit or the even cooler
   alternative to threefry (aesni). For the time being that's hardly
   worth it, threefry is already much better than the default RNG.
 */

#include "Random123/threefry.h"
#include "Random123/u01.h"

#define exa_rand(c,k) threefry2x32(c,k)

typedef threefry2x32_key_t randKey_t; 
typedef threefry2x32_ctr_t randCtr_t; 


class Randomness : public Checkpointable
{
public: 
  Randomness(randCtr_t seed); 
  randCtr_t generateSeed();     
  /** 
      @brief rebase a random number generator for a given generation 
   */ 
  void rebase(int num); 
  /** 
      @brief draw a floating point random  number in [0,1)
   */ 
  double drawRandDouble01(); 
  nat operator()() ; 
  double drawRandBiUnif(double x); 
  /** 
      @brief draw from a exponential distribution with parameter
      lambda
   */ 
  double drawRandExp(double lambda); 
  /**
     @brief draw a random number that is distributed uniformly on the
     exponential scale
   */ 
  double drawMultiplier(double multiplier); 
  /** 
      @brief draw a random number around param for a given window
   */ 
  double drawFromSlidingWindow(double param, double window); 
  
  randCtr_t getKey() const {return key; }
  /** 
      @brief draws a permutation 
   */ 
  void drawPermutation( int* perm, int n); 
  /**
     @brief draw r according to distribution given by weights. 
     
     NOTE sum of weights is not required to be 1.0
  */
  int drawSampleProportionally( double *weights, int numWeight ); 
  /** 
      @brief draws a random number form a dirichlet distribution
   */ 
  std::vector<double> drawRandDirichlet( const std::vector<double> &alphas); 
  /** 
      @brief draw a random number from a gamma distribution with shape
      parameter alpha and scale parameter beta
   */ 
  double drawRandGamma(double alpha, double beta); 
  /** 
      @brief This function should be called if the alphas for the dirichlet distribution are given
  */ 
  void drawRandDirichlet( std::vector<double> &results, const std::vector<double> &alphas); 
  /** 
      @brief This function should be called if the expected values for the dirichlet distribution are given      
  */
  void drawDirichletExpected(std::vector<double> &results, const std::vector<double> &mean,double scale); 
  /** 
      @brief draws integer uniformly from [0,n)
  */  
  int drawIntegerOpen(int upperBound); 
  /** 
      @brief draw integer uniformly from [0,n]
  */ 
  int drawIntegerClosed(int upperBound); 

  virtual void readFromCheckpoint( std::istream &in )  ; 
  virtual void writeToCheckpoint( std::ostream &out) const; 

  friend std::ostream& operator<<( std::ostream& out, const Randomness &rhs ); 

private:  		// METHODS
  void incrementNoLimit();

private:   			// ATTRIBUTES
  randCtr_t ctr; 
  randKey_t key; 
}; 

#endif
