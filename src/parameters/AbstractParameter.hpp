#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "ParameterContent.hpp"
#include "TreeAln.hpp"

class AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const = 0; 
  virtual ParameterContent extractParameter(const TreeAln &traln, std::vector<nat> model)  const  = 0;   

  void setSavedContent(ParameterContent& content) { savedContent = content; }
  ParameterContent& getSavedContent() {return savedContent; }
  
protected: 
  std::vector<nat> partitions; 
  ParameterContent savedContent; 
}; 

#endif
