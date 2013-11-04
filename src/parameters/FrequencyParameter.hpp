#ifndef FREQ_PARAMETER
#define FREQ_PARAMETER

#include "AbstractParameter.hpp"
#include "Category.hpp"
  
class FrequencyParameter : public AbstractParameter
{
public: 
  FrequencyParameter(nat id, nat idOfMyKind, std::vector<nat> partitions )
    : AbstractParameter(Category::FREQUENCIES, id , idOfMyKind, partitions)
  { 
  }

  FrequencyParameter(const FrequencyParameter &rhs )
    : AbstractParameter(rhs)
  {    
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const;   
  virtual AbstractParameter* clone () const {return new FrequencyParameter(*this); } 

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  ; 

  virtual void verifyContent(const TreeAln &traln, const ParameterContent &content) const;

  virtual void checkSanityPartitionsAndPrior(const TreeAln& traln) const ; 
}; 


#endif
