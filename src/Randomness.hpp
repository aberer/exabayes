#ifndef ___RANDOMNESS_H
#define ___RANDOMNESS_H

#include <iostream>
#include <vector>

#include "axml.h"

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

// class Branch; 

class Randomness
{
public: 
  Randomness(int seed);
  // ~Randomness(){assert(0); }
  randCtr_t generateSeed();   
  void rebase(int num){ ctr.v[1] = num; ctr.v[0] = 0; }
  
  int drawRandInt(int upperBound); 
  double drawRandDouble01(); 

  nat operator()() ; 

  double drawRandBiUnif(double x); 
  double drawRandExp(double lambda); 
  double drawMultiplier(double multiplier); 
  double drawFromSlidingWindow(double param, double window); 

  // branch drawSubtreeUniform(TreeAln &traln); 
  // branch drawInnerBranchUniform(TreeAln &traln); 

  void drawPermutation( int* perm, int n); 
  int drawSampleProportionally( double *weights, int numWeight ); 

  // Gamma(alpha, beta) sampling
  double drawRandGamma(double alpha, double beta); 

  //This function should be called if the alphas for the dirichlet distribution are given
  void drawRandDirichlet( std::vector<double> &results, const std::vector<double> &alphas); 


  //This function should be called if the expected values for the dirichlet distribution are given
  void drawDirichletExpected(std::vector<double> &results, const std::vector<double> &mean,double scale); 
  
  // Branch drawBranchUniform(TreeAln &traln); 
  // Branch drawInnerBranchUniform(const TreeAln &traln ); 
  int drawIntegerOpen(int upperBound); 

  /** @brief prints the RNG state  */ 
  friend std::ostream& operator<<( std::ostream& out, const Randomness &rhs ); 

  // Branch drawBranchWithInnerNode(const TreeAln &traln); 

  int drawIntegerClosed(int upperBound); 

  // nat drawInnerNode(const TreeAln &traln ); 

private:   
  randCtr_t ctr; 
  randKey_t key; 

  void incrementNoLimit();
  

  //density for dirichlet distribution with parameters "alphas" at point "values".

  // double gammaFunction(double alpha); 
  // double betaFunction(double *alpha, int length); 
  // void normalize(double * vector, int length, double normalizingConstant); 

}; 

#endif
