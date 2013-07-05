#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "ParameterContent.hpp"
#include "TreeAln.hpp"
#include "Priors.hpp"



enum class Category; 

class AbstractParameter
{
public: 
  
  AbstractParameter(Category cat, nat id)
    : id(id)
    , cat(cat) 
    // , modifiesBL (false)
  { }

  virtual AbstractParameter* clone () const =  0; 

  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const = 0; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  = 0;   

  void setSavedContent(const ParameterContent& content) { savedContent = content; }
  ParameterContent& getSavedContent() {return savedContent; }

  void addPartition(nat id){ partitions.push_back(id); }
  void setPrior(std::shared_ptr<AbstractPrior> _prior){prior = _prior; }

  Category getCategory() const {return cat; } 
  nat getId() const {return id; }
  std::vector<nat> getPartitions() const {return partitions; }
  AbstractPrior* getPrior() const { return prior.get(); }

  std::ostream&  printShort(std::ostream& out); 
  friend std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs); 
  
protected: 
  nat id; 
  Category cat; 
  std::vector<nat> partitions; 
  ParameterContent savedContent; 
  std::shared_ptr<AbstractPrior> prior; 
}; 

#endif