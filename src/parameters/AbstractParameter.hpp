#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "ParameterContent.hpp"
#include "TreeAln.hpp"
#include "priors/AbstractPrior.hpp"

enum class Category; 

class AbstractParameter
{
public: 
  
  AbstractParameter(Category cat, nat id)
    : id(id)
    , cat(cat) 
    , printToParamFile(true)
  {
  }

  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const = 0; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  = 0;   
  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln ) const = 0; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  = 0; 

  void addPartition(nat id){ partitions.push_back(id); }
  void setPrior(std::shared_ptr<AbstractPrior> _prior){prior = _prior; }

  Category getCategory() const {return cat; } 
  nat getId() const {return id; }
  std::vector<nat> getPartitions() const {return partitions; }
  AbstractPrior* getPrior() const { return prior.get(); }

  bool isPrintToParamFile() const {return printToParamFile; }

  std::ostream&  printShort(std::ostream& out); 
  friend std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs); 

  
  virtual AbstractParameter* clone() const = 0 ; 


protected: 
  nat id; 
  Category cat; 
  std::vector<nat> partitions; 
  std::shared_ptr<AbstractPrior> prior; 
  bool printToParamFile; 
}; 

#endif
