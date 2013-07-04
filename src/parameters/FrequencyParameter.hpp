#ifndef FREQ_PARAMETER
#define FREQ_PARAMETER

#include "AbstractParameter.hpp"
#include "Category.hpp"
  
class FrequencyParameter : public AbstractParameter
{
public: 
  FrequencyParameter(nat id)
    : AbstractParameter(Category::FREQUENCIES, id )
  {    
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const;   

}; 


#endif
