#ifndef _TOPOLOGY_PARAMETER
#define _TOPOLOGY_PARAMETER


#include "TreeAln.hpp"
#include "AbstractParameter.hpp"
#include "Category.hpp"

class TopologyParameter : public AbstractParameter
{
public: 
  TopologyParameter(nat id)
    : AbstractParameter(Category::TOPOLOGY, id)
  {
  }

  virtual void applyParameter(TreeAln& traln , const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new TopologyParameter(*this); } 
}; 

#endif
