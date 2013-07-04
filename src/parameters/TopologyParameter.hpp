#ifndef _TOPOLOGY_PARAMETER
#define _TOPOLOGY_PARAMETER


#include "TreeAln.hpp"
#include "AbstractParameter.hpp"

class TopologyParameter : public AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln, std::vector<nat> model)  const;   
  
}; 

#endif
