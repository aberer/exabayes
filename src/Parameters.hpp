#ifndef _PARAMETER_H
#define _PARAMETER_H

#include <cassert>

#include "GlobalVariables.hpp"
#include "TreeAln.hpp"
#include "PriorBelief.hpp"

class FrequencyParameter
{
public: 
  static void setParameters(TreeAln &traln, int model, vector<double> vals) 
  { 
    traln.setFrequenciesBounded(vals, model);
  }

  static vector<double> getParameters(TreeAln &traln, int model) { return traln.getFrequencies(model) ; }

  static void updatePrior(const TreeAln &traln, PriorBelief &prior, vector<double> oldVals, vector<double>newVals) 
  { 
    prior.updateFreq(oldVals, newVals); 
    prior.rescoreAllBranchLengths(traln); // TODO very INEFFICIENT ! 
  }

  static void init(TreeAln &traln, int model){ traln.initRevMat(model); }
  static category_t cat; 
}; 

class RevMatParameter
{
public: 
  static void setParameters(TreeAln &traln, int model, vector<double> values) 
  { 
    traln.setRevMatBounded(values, model) ; 
  } 
  
  static vector<double> getParameters(TreeAln &traln, int model) { return traln.getRevMat(model);  }  
  static void init(TreeAln &traln, int model){traln.initRevMat(model);}
  static void updatePrior(const TreeAln &traln, PriorBelief &prior, vector<double> oldVals, vector<double>newVals) 
  {  
    prior.updateRevMat(oldVals, newVals); 
    prior.rescoreAllBranchLengths(traln); // TODO verify INEFFICIENT ! 
  }
  
  static category_t cat; 
}; 


class RateHetParameter
{
public: 
  static void setParameters(TreeAln &traln, int model, vector<double> values) { assert(values.size() == 1) ; traln.setAlphaBounded(values[0],model);   }
  static vector<double> getParameters(TreeAln &traln, int model) { vector<double> result; result.push_back(traln.getAlpha(model)); return result;  }  
  static void init(TreeAln &traln, int model){traln.discretizeGamma(model); }
  static void updatePrior(const TreeAln &traln, PriorBelief &prior, vector<double> oldVals, vector<double>newVals) { prior.updateRateHet(oldVals[0], newVals[0]); }
  static category_t cat; 
}; 




#endif
