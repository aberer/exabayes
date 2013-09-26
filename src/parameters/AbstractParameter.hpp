#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "ParameterContent.hpp"

class AbstractPrior; 
class TreeAln; 

enum class Category; 

class AbstractParameter
{
public:   
  AbstractParameter(Category cat, nat id, nat idOfMyKind); 
  /** 
      @brief applies the parameter content to the tree 
   */ 
  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const = 0; 
  /** 
      @brief extracts the parameter 
   */ 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  = 0;   
  /** 
      @brief print a sample for this parameter 
   */ 
  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln ) const = 0; 
  /** 
      @brief print the names of all components of this parameter (e.g., the meaning of the various rates ) 
   */ 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  = 0; 
  /** 
      @brief adds a partition to the parameter (during setup)
   */ 
  void addPartition(nat id){ partitions.push_back(id); }
  /** 
      @brief sets the prior for this parameter 
   */ 
  void setPrior(std::shared_ptr<AbstractPrior> _prior){prior = _prior; }
  nat getIdOfMyKind() const {return idOfMyKind; }
  /** 
      @brief veriffies that content is compatible to this parameter (e.g., not too many rates). 

      This is a crude method merely for initialization (user input validation)
   */ 
  virtual void verifyContent(const TreeAln &traln, const ParameterContent &content) const  =  0; 

  ///////////////
  // OBSERVERS //
  ///////////////
  Category getCategory() const {return cat; } 
  nat getId() const {return id; }
  std::vector<nat> getPartitions() const {return partitions; }
  AbstractPrior* getPrior() const { return prior.get(); }
  bool isPrintToParamFile() const {return printToParamFile; }

  std::ostream&  printShort(std::ostream& out) const;  
  friend std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs); 
  virtual AbstractParameter* clone() const = 0 ; 

protected: 
  nat id; 
  nat  idOfMyKind;
  Category cat; 
  std::vector<nat> partitions; 
  std::shared_ptr<AbstractPrior> prior; 
  bool printToParamFile; 
}; 

#endif
