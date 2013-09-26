#include "AbstractProposal.hpp"


AbstractProposal::AbstractProposal( Category cat, const std::string& _name )  
  : name (_name)
  , category(cat)
  , needsFullTraversal(true)
  , inSetExecution(false)
{
} 


AbstractProposal::AbstractProposal( const AbstractProposal& rhs)
  : name(rhs.name)
  , sctr(rhs.sctr)
  , category(rhs.category)
  , relativeWeight(rhs.relativeWeight)
  , needsFullTraversal(rhs.needsFullTraversal)
  , inSetExecution(rhs.inSetExecution)
  , id(rhs.id)
{
  for(auto &v : rhs.primaryParameters)
    primaryParameters.emplace_back(v->clone()); 
  for(auto &v : rhs.secondaryParameters)
    secondaryParameters.emplace_back(v->clone()); 
}


void AbstractProposal::updateHastingsLog(double &hastings, double logValueToAdd, std::string whoDoneIt) 
{
  hastings += logValueToAdd; 
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


std::ostream&  operator<< ( std::ostream& out , const AbstractProposal& rhs) 
{
  out << rhs.name  << std::endl
      << "\tintegrating:\t"; 
  for(auto &r : rhs.primaryParameters)
    out << r.get() << ", "  ; 
  
  if(not rhs.secondaryParameters.empty() )
    {
      out << std::endl << "\talso modifying:\t" ; 
      for(auto &r : rhs.secondaryParameters ) 
	out << r.get() << ",\t" ; 
    }
  return out; 
}

 
void AbstractProposal::serialize( std::ostream &out)   const
{
  // auto && ss = std::stringstream{}; 
  // printShort(ss); 
  // std::string name = ss.str();
  // writeString(out, name); 

  cWrite(out, id);
  sctr.serialize(out) ; 
  writeToCheckpointCore(out); 
}


void AbstractProposal::deserialize( std::istream &in )
{
  // notice: name has already been read 
  sctr.deserialize(in); 
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


std::vector<nat> AbstractProposal::getAffectedPartitions() const 
{
  assert(primaryParameters.size() == 1); 
  return primaryParameters[0]->getPartitions();
} 
