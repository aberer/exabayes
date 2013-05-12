#include "PriorBelief.hpp"
#include "branch.h"




PriorBelief::PriorBelief()
  : logProb(0)
{
  
}


void PriorBelief::initPrior(const TreeAln& traln)
{
  logProb = scoreEverything(traln);
}


double PriorBelief::scoreEverything(const TreeAln &traln)
{
  double result = 0; 
  result += scoreBranchLengths(traln); 
  result += scoreRevMats(traln); 
  result += scoreRateHets(traln);
  result += scoreStateFreqs(traln);
  return result;
}



double PriorBelief::scoreBranchLengths(const TreeAln &traln)
{  
  double result = 0; 
  vector<branch> branches; 
  extractBranches(traln, branches);
  for(auto b : branches)
    {
      vector<double> tmp; 
      tmp.push_back(branchLengthToReal(traln.getTr(), b.length[0]));
      result += brPr->getLogProb(tmp);
    }
  return result; 
}


double PriorBelief::scoreRevMats(const TreeAln &traln)
{
  double result = 0; 
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      vector<double> rates; 
      for(int j = 0; j < 6; ++j )
	rates.push_back(traln.getSubstRate(i,j));
      result += revMatPr->getLogProb(rates);
    }
  return result; 
}

double PriorBelief::scoreRateHets(const TreeAln &traln)
{
  double result = 0; 
  for(int i = 0; i <  traln.getNumberOfPartitions(); ++i)
    {
      vector<double> alpha; 
      alpha.push_back(traln.getAlpha(i));
      result += rateHetPr->getLogProb(alpha);
    }
  return result; 
}


double PriorBelief::scoreStateFreqs(const TreeAln &traln)
{
  double result = 0; 
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      vector<double> freq; 
      for(int j = 0; j < 4; ++j)	// TODO 
	freq.push_back(traln.getFrequency(i,j));
      result += stateFreqPr->getLogProb(freq);
    }
  return result; 
}


void PriorBelief::verify(const TreeAln &traln)  
{
  assert(scoreEverything(traln) == logProb); 
}

