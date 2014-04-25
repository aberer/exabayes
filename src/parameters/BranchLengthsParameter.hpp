#ifndef _BRANCH_LENGTHS_PARAMETER  
#define _BRANCH_LENGTHS_PARAMETER  

#include "AbstractParameter.hpp"

#include "model/Category.hpp"

class BranchLengthsParameter : public AbstractParameter 
{
public: 
  
  BranchLengthsParameter(nat id, nat idOfMyKind, std::vector<nat> partitions)    
    : AbstractParameter(Category::BRANCH_LENGTHS, id, idOfMyKind, partitions,0 )
  {    
    _printToParamFile = false; 
  }

  BranchLengthsParameter(const BranchLengthsParameter& rhs)
    : AbstractParameter(rhs)
  {
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new BranchLengthsParameter(*this) ;  }

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const  {}
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  {} 

  virtual void verifyContent(const TreeAln&traln, const ParameterContent &content) const ; 

  virtual bool priorIsFitting(const AbstractPrior &prior, const TreeAln &traln) const; 

}; 


#endif
