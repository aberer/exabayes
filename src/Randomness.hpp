#ifndef ___RANDOMNESS_H
#define ___RANDOMNESS_H

#include "axml.h"
#include "branch.h"
#include <iostream>

typedef struct _branch  branch; 
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

#include <Random123/threefry.h>
#include <Random123/u01.h>

#define exa_rand(c,k) threefry2x32(c,k)

typedef threefry2x32_key_t randKey_t; 
typedef threefry2x32_ctr_t randCtr_t; 



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

  // branch drawSubtreeUniform(TreeAln &traln); 
  branch drawInnerBranchUniform(TreeAln &traln); 

  void drawPermutation( int* perm, int n); 
  int drawSampleProportionally( double *weights, int numWeight ); 

  // Gamma(alpha, beta) sampling
  double drawRandGamma(double alpha, double beta); 

  //This function should be called if the alphas for the dirichlet distribution are given
  void drawRandDirichlet( vector<double> &results, const vector<double> &alphas); 


  //This function should be called if the expected values for the dirichlet distribution are given
  void drawDirichletExpected(vector<double> &results, const vector<double> &mean,double scale); 

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
