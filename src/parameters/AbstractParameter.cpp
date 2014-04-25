#include "AbstractParameter.hpp" 
#include "system/extensions.hpp" 
#include "model/Category.hpp"


AbstractParameter::AbstractParameter(Category cat, nat id, nat idOfMyKind, std::vector<nat> partitions, nat paramPrio)
  : _id(id)
  , _idOfMyKind(idOfMyKind)
  , _cat(cat) 
  , _printToParamFile(true)
  , _partitions(partitions)
  , _paramPriority(paramPrio)
{
  assert(_partitions.size() > 0); 
}

AbstractParameter::AbstractParameter(const AbstractParameter& rhs)
  : _id(rhs._id)
  , _idOfMyKind(rhs._idOfMyKind)
  , _cat(rhs._cat)
  , _printToParamFile(rhs._printToParamFile)
  , _partitions(rhs._partitions)
  , _paramPriority(rhs._paramPriority)
{
  if(rhs._prior.get() != nullptr)
    _prior = std::unique_ptr<AbstractPrior>(rhs._prior->clone()); 
  assert(_partitions.size() > 0); 
}

std::ostream& operator<<(std::ostream &out, const AbstractParameter* rhs)
{
  return rhs->printShort(out);
}


std::ostream&  AbstractParameter::printShort(std::ostream& out) const
{
  out << CategoryFuns::getShortName(_cat) << "{" ; 

  formatRange(out, _partitions); 

  out << "}";     
  return out; 
}


void AbstractParameter::checkSanityPartitionsAndPrior(const TreeAln &traln) const 
{
}


void AbstractParameter::checkSanityPartitionsAndPrior_FreqRevMat(const TreeAln &traln) const 
{
  auto numStates = traln.getPartition(_partitions.at(0)).getStates();
  auto okay = bool{true}; 
  auto wrong = nat{0}; 
  for(auto p : _partitions)
    {
      if(numStates)
	okay &= traln.getPartition(p).getStates() == numStates; 
      if(not okay)
	{
	  wrong = p; 
	  break; 
	}
    }

  if(not okay)
    {
      std::cerr << "Error while processing parsed parameters: you tried to link " << _partitions[0] << " and " <<  wrong  << ". These partitions have a different number of states (probably DNA and PROT). Aborting." << std::endl; 
      exitFunction(-1, true); 
    }
}


bool AbstractParameter::priorIsFitting(const AbstractPrior &prior, const TreeAln &traln) const
{
  return true;   
}
