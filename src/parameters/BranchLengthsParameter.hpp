#ifndef _BRANCH_LENGTHS_PARAMETER  
#define _BRANCH_LENGTHS_PARAMETER  

#include "AbstractParameter.hpp"
#include "ComplexTuner.hpp"

#include "Category.hpp"

class BranchLengthsParameter : public AbstractParameter 
{
public: 			// inherited from SERIALIZABLE
  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const; 

public: 			// INHERITED METHODS 
  virtual void applyParameter(TreeAln& traln, const ParameterContent &content);
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new BranchLengthsParameter(*this) ;  }

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const  {}
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  {} 

  virtual void verifyContent(const TreeAln&traln, const ParameterContent &content) const ; 

  virtual bool priorIsFitting(const AbstractPrior &prior, const TreeAln &traln) const; 

  virtual ParamAttribute getAttributes() const {  return { _convTuner , _nonConvTuner} ;  } 
  virtual void setAttributes(ParamAttribute attr) { _convTuner = attr._convTuner;  _nonConvTuner = attr._nonConvTuner;  }
  
  virtual double getMeanSubstitutionRate() const ;
  virtual void updateMeanSubstRate(const TreeAln& traln) ;  
  virtual void setMeanSubstitutionRate(double fac) {_fracChange = fac; }

  virtual log_double getPriorValue(const TreeAln& traln) const { assert(0);  return log_double::fromAbs(1);  }

public: 			// METHODS 
  BranchLengthsParameter(nat id, nat idOfMyKind, std::vector<nat> partitions)    ; 
  BranchLengthsParameter(const BranchLengthsParameter& rhs) = default; 

private: 			// ATTRIBUTES 
  ComplexTuner _convTuner; // convergence parameter for distribution proposals  
  ComplexTuner _nonConvTuner; // non-convergence parameter for the distribution proposals 
  double _fracChange; 
}; 


#endif
