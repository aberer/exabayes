#ifndef FREQ_PARAMETER
#define FREQ_PARAMETER

#include "AbstractParameter.hpp"
  
class FrequencyParameter : public AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln, std::vector<nat> model)  const;   

}; 


#endif
