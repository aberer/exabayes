#ifndef RATE_HET_PARAMETER
#define RATE_HET_PARAMETER

#include "AbstractParameter.hpp"
#include "Category.hpp"
  
class RateHetParameter : public AbstractParameter
{
public: 

  RateHetParameter(nat id, nat idOfMyKind )
    : AbstractParameter(Category::RATE_HETEROGENEITY, id, idOfMyKind )
  {
  }
  
  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new RateHetParameter(*this); } 

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const ; 

  virtual void verifyContent(const TreeAln&traln, const ParameterContent &content) const; 
}; 

#endif
