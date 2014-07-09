#include "AbstractProposal.hpp"


AbstractProposal::AbstractProposal( Category cat, std::string name,double weight, double minTuning, double maxTuning, bool needsFullTraversal )  
  : _name(name)
  , _category(cat)
  , _relativeWeight(weight)
  , _needsFullTraversal(needsFullTraversal)
  , _inSetExecution(false)
  , _forProteinOnly(false)
  , _minTuning{minTuning}
  , _maxTuning{maxTuning}
  , _numTaxNeeded{4}
  , _usingOptimizedBranches(false)
{
} 


AbstractProposal::AbstractProposal( const AbstractProposal& rhs)
  : _sctr(rhs._sctr)
  , _category(rhs._category)
  , _relativeWeight(rhs._relativeWeight)
  , _needsFullTraversal(rhs._needsFullTraversal)
  , _inSetExecution(rhs._inSetExecution)
  , _id(rhs._id)
  , _forProteinOnly(rhs._forProteinOnly)
  , _minTuning{rhs._minTuning}
  , _maxTuning{rhs._maxTuning}
  , _numTaxNeeded{rhs._numTaxNeeded}
  , _usingOptimizedBranches(rhs._usingOptimizedBranches)
{
  this->_name = rhs._name; 

  for(auto &v : rhs._primaryParameters)
    _primaryParameters.emplace_back(v->clone()); 
  for(auto &v : rhs._secondaryParameters)
    _secondaryParameters.emplace_back(v->clone()); 
}


// void AbstractProposal::updateHastingsLog(double &hastings, double logValueToAdd, std::string whoDoneIt) 
// {
//   hastings += logValueToAdd; 
// }


std::ostream& AbstractProposal::printShort(std::ostream &out)  const 
{
  out << _name << "( " ;  
    
  bool isFirst = true; 
  for(auto &v : _primaryParameters)
    {
      if(not isFirst)
	out << ","; 
      else 
	isFirst = false; 
      v->printShort(out); 
    }

  if(_secondaryParameters.size() > 0)
    {
      out << ";"; 
      isFirst = true; 
      for(auto &v : _secondaryParameters)
	{
	  if(not isFirst)
	    out << ","; 
	  else 
	    isFirst = false; 
	  v->printShort(out); 
	}
    }

  printParams(out);

  out << " )"; 
  return out; 
}


std::ostream& AbstractProposal::printNamePartitions(std::ostream &out)
{
  out << _name  << "(" ; 
  assert(_primaryParameters.size() == 1); 
  bool isFirst= true; 
  for (auto v : _primaryParameters[0]->getPartitions()) 
    {
      if( not isFirst)
	out << ","; 
      else 
	isFirst = false; 
      out << v ; 
    }
  out << ")" ; 
  return out; 
}


std::ostream&  operator<< ( std::ostream& out , const AbstractProposal& rhs) 
{
  out << rhs._name  << std::endl
      << "\tintegrating:\t"; 
  for(auto &r : rhs._primaryParameters)
    out << r.get() << ", "  ; 
  
  if(not rhs._secondaryParameters.empty() )
    {
      out << std::endl << "\talso modifying:\t" ; 
      for(auto &r : rhs._secondaryParameters ) 
	out << r.get() << ",\t" ; 
    }
  return out; 
}


void AbstractProposal::serialize( std::ostream &out)   const
{
  _sctr.serialize(out) ; 
  writeToCheckpointCore(out); 
}


void AbstractProposal::deserialize( std::istream &in )
{
  // auto id = cRead<decltype(_id)>(in); 
  // std::cout << "read id " << id << ", expected " << _id << std::endl; 
  // assert(id == _id); 
  _sctr.deserialize(in); 
  readFromCheckpointCore(in); 
}


std::vector<AbstractParameter*> AbstractProposal::getPrimaryParameterView() const
{
  auto result = std::vector<AbstractParameter*> {}; 
  for(auto &v : _primaryParameters)
    result.push_back(v.get()); 
  return result; 
}

 
std::vector<AbstractParameter*> AbstractProposal::getSecondaryParameterView() const 
{
  std::vector<AbstractParameter*> result; 
  for(auto &v : _secondaryParameters)
    result.push_back(v.get()); 
  return result; 
}


std::vector<AbstractParameter*> AbstractProposal::getBranchLengthsParameterView() const 
{
  auto result = std::vector<AbstractParameter*> {}; 
  for(auto &p : _primaryParameters)
    if(p->getCategory() == Category::BRANCH_LENGTHS)
      result.push_back(p.get()); 

  for(auto &p : _secondaryParameters ) 
    if(p->getCategory() == Category::BRANCH_LENGTHS)
      result.push_back(p.get());
  return result; 
} 


std::vector<nat> AbstractProposal::getAffectedPartitions() const 
{
  assert(_primaryParameters.size() == 1); 
  return _primaryParameters[0]->getPartitions();
} 


std::array<bool,3>
AbstractProposal::getBranchProposalMode() const 
{
  bool outer = false; 
  bool multiply = false; 
  bool sequential = false; 

  auto branchOpt = std::getenv("PROPOSE_BRANCHES"); 
  if(branchOpt != NULL )
    {
      auto &&iss =  std::istringstream(branchOpt); 
      int mode = 0; 
      iss >> mode; 
      if(mode == 0)
	{
	}
      else if(mode == 1)
	{
	  multiply = true ; 
	}
      else if(mode == 2)
	{
	  multiply = true; 
	  outer = true; 
	}
      else if(mode == 3 )
	{
	  multiply = true; 
	  sequential = true; 
	}
      else if(mode == 4 )
	{
	  multiply = true; 
	  outer = true; 
	  sequential = true; 
	}
      else 
	assert(0);
    }

  return {{multiply, outer, sequential }}; 
}


/**
   @brief log tuning for a parameter 
   
   @return tuned parameter    
 */
double AbstractProposal::tuneParameter(int batch, double accRatio, double parameter, bool inverse )
{  
  double delta = 1.0 / sqrt(batch + 1 );
  // delta = 0.1 < delta ? 0.1 : delta;

  assert(_minTuning != 0  || _maxTuning != 0 ); 

  double logTuning = log(parameter);
  
  if(inverse)
    logTuning += (accRatio > TARGET_RATIO)  ? -delta : +delta ;
  else 
    logTuning += (accRatio > TARGET_RATIO)  ? +delta : -delta ;
  
  double newTuning = exp(logTuning);

  if (_minTuning <  newTuning && newTuning < _maxTuning)
    return newTuning; 
  else 
    return parameter; 
}
