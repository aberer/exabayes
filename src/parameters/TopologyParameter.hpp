#ifndef _TOPOLOGY_PARAMETER
#define _TOPOLOGY_PARAMETER


#include "TreeAln.hpp"
#include "AbstractParameter.hpp"
#include "Category.hpp"

class TopologyParameter : public AbstractParameter
{
public: 
  TopologyParameter(nat id, nat idOfMyKind )
    : AbstractParameter(Category::TOPOLOGY, id, idOfMyKind)
  {
    printToParamFile = false; 
  }

  virtual void applyParameter(TreeAln& traln , const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new TopologyParameter(*this); } 

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const {}
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  {} 
}; 

#endif
