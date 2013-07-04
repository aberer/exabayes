#ifndef RATE_HET_PARAMETER
#define RATE_HET_PARAMETER

#include "AbstractParameter.hpp"
  
class RateHetParameter : public AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln, std::vector<nat> model)  const;   

}; 

#endif
