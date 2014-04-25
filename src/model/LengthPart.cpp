#include "LengthPart.hpp"
#include "TreeAln.hpp"
#include "parameters/AbstractParameter.hpp"


void LengthPart<double>::extractLength(const TreeAln &traln, const BranchPlain& branch, const AbstractParameter*  param)
{
  auto p = branch.findNodePtr(traln); 
  length = p->z[param->getPartitions()[0]]; 
}


void LengthPart<std::vector<double>>::extractLength(const TreeAln &traln, const BranchPlain &branch, const std::vector<AbstractParameter*> &params)
{
  auto p = branch.findNodePtr(traln); 
  lengths.clear();
  for(auto &param: params )
    lengths.push_back(p->z[param->getPartitions()[0]]) ; 
}


double LengthPart<double>::getInterpretedLength(const TreeAln &traln, const AbstractParameter* param) const
{ 
  double fracC = traln.getMeanSubstitutionRate(param->getPartitions()); 
  return -log(length) * fracC; 
}


void LengthPart<double>::setConvertedInternalLength(const TreeAln& traln, const AbstractParameter* param, double length) 
{
  double fracC = traln.getMeanSubstitutionRate(param->getPartitions());  
  double internalLength = exp(- length / fracC); 
  this->length = internalLength; 
} 




double LengthPart<std::vector<double> >::getLength (const AbstractParameter* param) const 
{
  return lengths.at(param->getIdOfMyKind()) ;
}



