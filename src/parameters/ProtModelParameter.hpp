#ifndef PROT_MODEL_PARAMETER
#define PROT_MODEL_PARAMETER

#include "model/Category.hpp"
#include "parameters/AbstractParameter.hpp"

class ProtModelParameter : public AbstractParameter
{
public: 
  ProtModelParameter(nat id, nat idOfMyKind, std::vector<nat> partitions)
    : AbstractParameter(Category::AA_MODEL, id, idOfMyKind, partitions,1)
  {
  }
  
  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const ; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  ;   
  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln ) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  ; 
  virtual void verifyContent(const TreeAln &traln, const ParameterContent &content) const  ; 
  virtual AbstractParameter* clone() const { return new ProtModelParameter(*this); } 

  virtual void checkSanityPartitionsAndPrior(const TreeAln &traln) const  ;

private: 
  
  
}; 


#endif
