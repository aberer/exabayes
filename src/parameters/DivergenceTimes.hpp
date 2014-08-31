#ifndef DIVERGENCETIME_H
#define DIVERGENCETIME_H


#include "AbstractParameter.hpp"

#include "NodeAge.hpp"

class DivergenceTimes : public AbstractParameter
{
public:
  DivergenceTimes(nat id, nat idOfMyKind, std::vector<nat> partitions, nat numberOfTaxa);
  virtual ~DivergenceTimes(){}

  virtual void applyParameter(TreeAln& traln,  const ParameterContent &content) const ; 
  virtual ParameterContent extractParameter(const TreeAln &traln)  const  ;   
  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln ) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  ; 
  virtual void verifyContent(const TreeAln &traln, const ParameterContent &content) const; 

  
  virtual AbstractParameter* clone() const  
  {
    return new DivergenceTimes(*this); 
  } 


private: 
  std::vector<NodeAge> _nodeAges; 

};



#endif /* DIVERGENCETIME_H */
