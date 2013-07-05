#ifndef RATE_HET_PARAMETER
#define RATE_HET_PARAMETER

#include "AbstractParameter.hpp"
#include "Category.hpp"
  
class RateHetParameter : public AbstractParameter
{
public: 

  RateHetParameter(nat id)
    : AbstractParameter(Category::RATE_HETEROGENEITY, id)
  {
    // modifiesBL = false; 
  }
  
  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new RateHetParameter(*this); } 
}; 

#endif
