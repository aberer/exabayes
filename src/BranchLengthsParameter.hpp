#ifndef _BRANCH_LENGTHS_PARAMETER  
#define _BRANCH_LENGTHS_PARAMETER  

class BranchLengthsParameter : public AbstractParameter 
{
public: 
  virtual void applyParameter(TreeAln& traln, vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent& extractParameter(const TreeAln &traln, vector<nat> model)  const;   

}; 


#endif
