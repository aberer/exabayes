#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "ParameterContent.hpp"

class AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, vector<nat> model, ParameterContent &content) const = 0; 
  virtual ParameterContent& extractParameter(const TreeAln &traln, vector<nat> model)  const  = 0;   

  void setSavedContent(ParameterContent& content) { savedContent = content; }
  ParameterContent& getSavedContent() {return savedContent; }
  
protected: 
  vector<nat> partitions; 
  ParameterContent savedContent; 
}; 

#endif
