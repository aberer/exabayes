#ifndef DIVERGENCERATE_H
#define DIVERGENCERATE_H

#include "AbstractParameter.hpp"



class DivergenceRates : public AbstractParameter
{
public:
  DivergenceRates(nat id, nat idOfMyKind, std::vector<nat> partitions, nat numberOfTaxa, nat numRateCats); 
  virtual ~DivergenceRates(){}

  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const ; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  ;   
  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln ) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  ; 
  virtual void verifyContent(const TreeAln &traln, const ParameterContent &content) const ; 
  
  virtual AbstractParameter* clone() const  
  {
    return new DivergenceRates(*this); 
  } 

private: 
  std::vector<int> _rateAssignments; 
  std::vector<double> _rates; 
};




#endif /* DIVERGENCERATE_H */
