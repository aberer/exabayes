#ifndef REV_MAT_PARAMETER
#define REV_MAT_PARAMETER


#include "Category.hpp"
#include "AbstractParameter.hpp"
  
class RevMatParameter : public AbstractParameter
{
public: 
  RevMatParameter(nat id) 
    : AbstractParameter(Category::SUBSTITUTION_RATES, id)
  {
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   

}; 

#endif
