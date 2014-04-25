#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "model/TreeAln.hpp"
#include "ParameterContent.hpp"
#include "priors/AbstractPrior.hpp"

enum class Category; 

class AbstractParameter
{
public:   
  AbstractParameter(Category cat, nat id, nat idOfMyKind, std::vector<nat> partitions, nat paramPrio); 
  AbstractParameter(const AbstractParameter& rhs); 

  /** 
      @brief applies the parameter content to the tree 
   */ 
  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const = 0; 
  virtual void applyParameterRaw(TreeAln &traln, const ParameterContent & content) const {}
  /** 
      @brief extracts the parameter 
   */ 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  = 0;   
  virtual ParameterContent extractParameterRaw(const TreeAln& traln) const {return ParameterContent{}; }
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
  void addPartition(nat id){ _partitions.push_back(id); }
  /** 
      @brief sets the prior for this parameter 
   */ 
  void setPrior(const std::unique_ptr<AbstractPrior> &prior){_prior = std::unique_ptr<AbstractPrior>(prior->clone()); }
  nat getIdOfMyKind() const {return _idOfMyKind; }
  /** 
      @brief veriffies that content is compatible to this parameter (e.g., not too many rates). 

      This is a crude method merely for initialization (user input validation)
   */ 
  virtual void verifyContent(const TreeAln &traln, const ParameterContent &content) const  =  0; 

  ///////////////
  // OBSERVERS //
  ///////////////
  Category getCategory() const {return _cat; } 
  nat getId() const {return _id; }
  std::vector<nat> getPartitions() const {return _partitions; }
  AbstractPrior* getPrior() const { return _prior.get(); }
  bool isPrintToParamFile() const {return _printToParamFile; }

  std::ostream&  printShort(std::ostream& out) const;  
  friend std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs); 
  virtual AbstractParameter* clone() const = 0 ; 


  virtual bool priorIsFitting(const AbstractPrior &prior, const TreeAln &traln) const; 

  virtual void checkSanityPartitionsAndPrior(const TreeAln &traln) const ; 

  nat getParamPriority() const {return _paramPriority; } 

protected: 			// METHODS
  void checkSanityPartitionsAndPrior_FreqRevMat(const TreeAln &traln) const ; 

protected: 			// ATTRIBUTES
  nat _id; 
  nat _idOfMyKind;
  Category _cat; 
  std::unique_ptr<AbstractPrior> _prior; 
  bool _printToParamFile; 
  std::vector<nat> _partitions; 
  nat _paramPriority; 
}; 

#endif
