#include "AbstractProposal.hpp"


AbstractProposal::AbstractProposal( const AbstractProposal& rhs)
  : name(rhs.name)
  , sctr(rhs.sctr)
  , category(rhs.category)
  , relativeWeight(rhs.relativeWeight)
  , needsFullTraversal(rhs.needsFullTraversal)
  , inSetExecution(rhs.inSetExecution)
{
  for(auto &v : rhs.primaryParameters)
    primaryParameters.emplace_back(v->clone()); 
  for(auto &v : rhs.secondaryParameters)
    secondaryParameters.emplace_back(v->clone()); 
}


/**
   @brief valToAdd must not be on the log-scale 
 */ 
void AbstractProposal::updateHastings(double &hastings, double valToAdd, std::string whoDoneIt) 
{
#ifdef DEBUG_HASTINGS  
  if(whoDoneIt.compare("branchCollapser" ) == 0 )
    tout << setprecision(6) << whoDoneIt << " updates hastings " << hastings << " with " << valToAdd ; 
#endif

  hastings += log(valToAdd); 	// we are logarithmic now   

#ifdef DEBUG_HASTINGS
  if(whoDoneIt.compare("branchCollapser") == 0)
    tout <<  " => " << hastings << endl; 
#endif
}


std::ostream& AbstractProposal::printShort(std::ostream &out)  const 
{
  out << this->name << "( " ;  
    
  bool isFirst = true; 
  for(auto &v : primaryParameters)
    {
      if(not isFirst)
	out << ","; 
      else 
	isFirst = false; 
      v->printShort(out); 
    }

  if(secondaryParameters.size() > 0)
    {
      out << ";"; 
      isFirst = true; 
      for(auto &v : secondaryParameters)
	{
	  if(not isFirst)
	    out << ","; 
	  else 
	    isFirst = false; 
	  v->printShort(out); 
	}
    }
  out << " )"; 
  return out; 
}


std::ostream& AbstractProposal::printNamePartitions(std::ostream &out)
{
  out << name  << "(" ; 
  assert(primaryParameters.size() == 1); 
  bool isFirst= true; 
  for (auto v : primaryParameters[0]->getPartitions()) 
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


std::ostream&  operator<< ( std::ostream& out , AbstractProposal* rhs) 
{
  out << rhs->name  << std::endl
      << "\tintegrating:\t"; 
  for(auto &r : rhs->primaryParameters)
    out << r.get() << ", "  ; 
  
  if(not rhs->secondaryParameters.empty() )
    {
      out << std::endl << "\talso modifying:\t" ; 
      for(auto &r : rhs->secondaryParameters ) 
	out << r.get() << ",\t" ; 
    }
  return out; 
}

 
void AbstractProposal::writeToCheckpoint( std::ostream &out)   const
{
  std::stringstream ss ; 
  printShort(ss); 
  std::string name = ss.str();
  writeString(out, name); 
  sctr.writeToCheckpoint(out) ; 
  writeToCheckpointCore(out); 
}


void AbstractProposal::readFromCheckpoint( std::istream &in )
{
  // notice: name has already been read 
  sctr.readFromCheckpoint(in); 
  readFromCheckpointCore(in); 
}


std::vector<AbstractParameter*> AbstractProposal::getPrimaryParameterView() const
{
  std::vector<AbstractParameter*> result; 
  for(auto &v : primaryParameters)
    result.push_back(v.get()); 
  return result; 
}

 
std::vector<AbstractParameter*> AbstractProposal::getSecondaryParameterView() const 
{
  std::vector<AbstractParameter*> result; 
  for(auto &v : secondaryParameters)
    result.push_back(v.get()); 
  return result; 
}


std::vector<AbstractParameter*> AbstractProposal::getBranchLengthsParameterView() const 
{
  auto result = std::vector<AbstractParameter*> {}; 
  for(auto &p : primaryParameters)
    if(p->getCategory() == Category::BRANCH_LENGTHS)
      result.push_back(p.get()); 

  for(auto &p : secondaryParameters ) 
    if(p->getCategory() == Category::BRANCH_LENGTHS)
      result.push_back(p.get());
  return result; 
} 
