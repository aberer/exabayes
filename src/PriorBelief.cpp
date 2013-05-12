#include "PriorBelief.hpp"
#include "branch.h"
#include "GlobalVariables.hpp"
#include "Priors.hpp"




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
  double verified = scoreEverything(traln); 

  if(fabs (verified - logProb)  > ACCEPTED_LIKELIHOOD_EPS) 
    {
      cerr << "WARNING: lnPrior was " <<  logProb << " while verification says, it should be " << verified << "\t"
	   << "DIFF:\t" << verified - logProb  << endl; 
      assert(0);       
    }  
}

// TODO some defaults may change in the light of various run
// configurations => we'd need a run-configuration object for that
// this would also solve our problem with various states and so on
void PriorBelief::addStandardPriors()
{
  if( brPr.get() == nullptr)
    {
      setBranchLengthPrior(shared_ptr<ExponentialPrior>(new ExponentialPrior(10)));
      tout << "No branch length prior specified. Falling back to default " << brPr.get() << " as branch length prior." << endl; 
    }

  if( revMatPr.get() == nullptr)
    {
      vector<double> tmp ; 
      for(int i = 0; i < 6; ++i)
	tmp.push_back(1.0); 
      setRevMatPrior(shared_ptr<DirichletPrior>(new DirichletPrior(tmp ))); 
      tout << "No prior for the reversible substitution matrix specified. Falling back to default "<<  revMatPr.get() << " as prior for the reversible matrix." << endl; 
    }

  if(rateHetPr.get() == nullptr)
    {
      setRateHetPrior(shared_ptr<UniformPrior>(new UniformPrior(1e-6, 200))); 
      tout << "No prior for rate heterogeneity specified. Falling back to default "<< rateHetPr.get() <<  " as prior for rate heterogeneity among sites." << endl; 
    }

  if(stateFreqPr.get() == nullptr)
    {
      vector<double> tmp ; 
      for(int i = 0; i < 4; ++i)
	tmp.push_back(1.0); 
      setStateFreqPrior(shared_ptr<DirichletPrior>(new DirichletPrior(tmp))); 
      tout << "No prior for state frequencies specified. Falling back to default " << stateFreqPr.get() << " as prior for state frequencies. " << endl; 
    }  
}



void PriorBelief::updateBranchLength(double oldValue, double newValue)
{
  vector<double> oldV = {oldValue}; 
  vector<double> newV = {newValue};   
  logProb +=  ( brPr->getLogProb(newV) - brPr->getLogProb(oldV) ) ; 
} 

void PriorBelief::updateRevMat( vector<double> oldValues, vector<double> newValues)
{
  logProb += revMatPr->getLogProb(newValues) - revMatPr->getLogProb(oldValues); 
}
 
void PriorBelief::updateFreq(vector<double> oldValues, vector<double> newValues)
{
  logProb += stateFreqPr->getLogProb(newValues) - revMatPr->getLogProb(oldValues); 
}
 
void PriorBelief::updateRateHet(double oldValue, double newValue)
{
  vector<double> oldV = {oldValue}; 
  vector<double> newV = {newValue} ;
  logProb += (rateHetPr->getLogProb(newV) - rateHetPr->getLogProb(oldV)); 
} 
