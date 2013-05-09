#ifndef _PARAMETER_H
#define _PARAMETER_H

#include "TreeAln.hpp"
#include "categoryType.h"




class FrequencyParameter
{
public: 
  static void setParameters(TreeAln &traln, int model, vector<double> vals) {     for(nat i = 0; i < vals.size(); ++i) traln.setFrequencySave(vals[i], model, i); }
  static vector<double> getParameters(TreeAln &traln, int model) { vector<double> result; for(int i = 0; i < 4; ++i) result.push_back( traln.getFrequency(model, i) ); return result; }
  static void init(TreeAln &traln, int model){
    // assert(0);
    traln.initRevMat(model);
  }

  static category_t cat; 
}; 


class RevMatParameter
{
public: 
  static void   setParameters(TreeAln &traln, int model, vector<double> values) { for(nat i = 0; i < values.size() ;++i) traln.setSubstSave(values[i], model, i);   }
  static vector<double> getParameters(TreeAln &traln, int model) { vector<double> result; for(int i = 0; i < 6 ; ++i) result.push_back(traln.getSubstRate(model, i)); return result;  }  
  static void init(TreeAln &traln, int model){traln.initRevMat(model);}
  static category_t cat; 
}; 


// nasty: position of course is not used 
class RateHetParameter
{
public: 
  static void setParameters(TreeAln &traln, int model, vector<double> values) { assert(values.size() == 1) ; traln.setAlphaSave(values[0],model);   }
  static vector<double> getParameters(TreeAln &traln, int model) { vector<double> result; result.push_back(traln.getAlpha(model)); return result;  }  
  static void init(TreeAln &traln, int model){traln.discretizeGamma(model); }
  static category_t cat; 
}; 




#endif
