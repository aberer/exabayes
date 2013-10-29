#include "DiscreteModelPrior.hpp"
#include "ProtModel.hpp"
// #include "common.h"


DiscreteModelPrior:: DiscreteModelPrior(std::unordered_map<ProtModel,double> model)
  :_modelProbs{model}
{
} 


ParameterContent DiscreteModelPrior::getInitialValue() const 
{
  auto result = ParameterContent{}; 
  auto aModel = std::get<0>( *(_modelProbs.begin() )); 
  result.protModel.push_back( aModel);  
  return result;
} 


double DiscreteModelPrior::accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param , double myOld, double myNew ) const 
{
  assert(0); 
}

 
std::vector<double> DiscreteModelPrior::drawFromPrior(Randomness &rand)  const 
{
  assert(0);
  // do we still need this routine? 
}

double DiscreteModelPrior::getLogProb( const ParameterContent& content ) const 
{
  assert(content.values.size() == 1 ); 
  
  auto model = ProtModel(content.values[0]); 
  assert(_modelProbs.find(model) != _modelProbs.end()); 

  return log(_modelProbs.at(model));
}

 
void DiscreteModelPrior::print(std::ostream &out) const 
{
  if(_modelProbs.size() == 1 )
    {
      auto model = std::get<0>(* _modelProbs.begin()); 
      out << "Fixed(" << ProtModelFun::getName(model) << ")" ; 
    }
  else 
    {
      out << "Discrete(" ; 
      auto isFirst = bool{true}; 
      for(auto &modelPair : _modelProbs)
	{
	  if(isFirst)
	    isFirst = false; 
	  else 
	    out << "," ;

	  out << ProtModelFun::getName(std::get<0>(modelPair)) << "=" << SOME_FIXED_PRECISION <<  std::get<1>(modelPair) ; 
	}
      out << ")"; 
    }
}

