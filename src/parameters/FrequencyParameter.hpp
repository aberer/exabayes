#ifndef FREQ_PARAMETER
#define FREQ_PARAMETER

#include "AbstractParameter.hpp"
#include "Category.hpp"
  
class FrequencyParameter : public AbstractParameter
{
public: 
  FrequencyParameter(nat id )
    : AbstractParameter(Category::FREQUENCIES, id )
  { 
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const;   
  virtual AbstractParameter* clone () const {return new FrequencyParameter(*this); } 

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  ; 
}; 


#endif
