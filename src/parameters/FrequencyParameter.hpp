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
    // modifiesBL = true; 
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const;   
  virtual AbstractParameter* clone () const {return new FrequencyParameter(*this); } 

}; 


#endif
