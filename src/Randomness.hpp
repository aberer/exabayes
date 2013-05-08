#ifndef ___RANDOMNESS_H
#define ___RANDOMNESS_H

#include "rng.h"
#include "axml.h"
#include "branch.h"
#include <iostream>

typedef struct _branch  branch; 
class TreeAln; 



class Randomness
{
public: 
  Randomness(int seed);
  randCtr_t generateSeed();   
  void rebase(int num){ ctr.v[1] = num; ctr.v[0] = 0; }

  int drawRandInt(int upperBound); 
  double drawRandDouble01(); 

  double drawRandBiUnif(double x); 
  double drawRandExp(double lambda); 
  double drawMultiplier(double multiplier); 
  double drawFromSlidingWindow(double param, double window); 

  branch drawSubtreeUniform(TreeAln &traln); 
  branch drawInnerBranchUniform(TreeAln &traln); 

  void drawPermutation( int* perm, int n); 
  int drawSampleProportionally( double *weights, int numWeight ); 

  // Gamma(alpha, beta) sampling
  double drawRandGamma(double alpha, double beta); 

  //This function should be called if the alphas for the dirichlet distribution are given
  void drawRandDirichlet( double* results, double* alphas,  int length); 


  //This function should be called if the expected values for the dirichlet distribution are given
  void drawDirichletExpected(double* results, double* mean,double scale, int length); 

  branch drawBranchUniform(TreeAln &traln); 

  /** @brief prints the RNG state  */ 
  friend ostream& operator<<( ostream& out, const Randomness &rhs ); 

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
