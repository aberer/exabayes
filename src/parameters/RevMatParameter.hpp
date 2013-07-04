#ifndef REV_MAT_PARAMETER
#define REV_MAT_PARAMETER


#include "AbstractParameter.hpp"
  
class RevMatParameter : public AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln, std::vector<nat> model)  const;   

}; 

#endif
