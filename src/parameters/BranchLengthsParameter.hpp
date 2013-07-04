#ifndef _BRANCH_LENGTHS_PARAMETER  
#define _BRANCH_LENGTHS_PARAMETER  

#include "AbstractParameter.hpp"

class BranchLengthsParameter : public AbstractParameter 
{
public: 
  virtual void applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln, std::vector<nat> model)  const;   
}; 


#endif
