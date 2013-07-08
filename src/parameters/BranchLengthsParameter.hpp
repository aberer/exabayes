#ifndef _BRANCH_LENGTHS_PARAMETER  
#define _BRANCH_LENGTHS_PARAMETER  

#include "AbstractParameter.hpp"

#include "Category.hpp"

class BranchLengthsParameter : public AbstractParameter 
{
public: 
  
  BranchLengthsParameter(nat id )    
    : AbstractParameter(Category::BRANCH_LENGTHS, id )
  {    
    printToParamFile = false; 
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new BranchLengthsParameter(*this) ;  }


  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const  {}
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  {} 
}; 


#endif
