#include "LengthPart.hpp"
#include "TreeAln.hpp"
#include "AbstractParameter.hpp"


bool LengthPart<double>::lenEqual(const LengthPart<double> &rhs) const 
{
  return std::fabs(length  - rhs.length) < std::numeric_limits<double>::epsilon();
}

bool LengthPart<std::vector<double> >::lenEqual(const LengthPart<std::vector<double> > &rhs) const 
{
  auto result = true; 
  result &= rhs.lengths.size() == lengths.size(); 
  for(auto i  = 0u; i < rhs.lengths.size() ; ++i)
    result &= std::fabs(rhs.lengths[i] - lengths[i] ) < std::numeric_limits<double>::epsilon();
  
  return result; 
}



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


double LengthPart<double>::getInterpretedLength(const AbstractParameter* param) const
{ 
  double fracC = param->getMeanSubstitutionRate();
  return -log(length) * fracC; 
}


// meh, could better go to branch length parameter 
void LengthPart<double>::setConvertedInternalLength(const AbstractParameter* param, double len) 
{
  double fracC = param->getMeanSubstitutionRate();
  double internalLength = exp(- len / fracC); 
  this->length = internalLength; 
} 




double LengthPart<std::vector<double> >::getLength (const AbstractParameter* param) const 
{
  return lengths.at(param->getIdOfMyKind()) ;
}



