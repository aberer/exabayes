#include "AbstractProposal.hpp"


std::vector<AbstractParameter*> AbstractProposal::getPrimVar() const
{
  std::vector<AbstractParameter*> result; 
  for(auto &v : primVar)
    result.push_back(v.get()); 
  return result; 
}

 
std::vector<AbstractParameter*> AbstractProposal::getSecVar() const 
{
  std::vector<AbstractParameter*> result; 
  for(auto &v : secVar)
    result.push_back(v.get()); 
  return result; 
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
  for(auto &v : primVar)
    {
      if(not isFirst)
	out << ","; 
      else 
	isFirst = false; 
      v->printShort(out); 
    }

  if(secVar.size() > 0)
    {
      out << ";"; 
      isFirst = true; 
      for(auto &v : secVar)
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
  assert(primVar.size() == 1); 
  bool isFirst= true; 
  for (auto v : primVar[0]->getPartitions()) 
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


std::ostream&  operator<< ( std::ostream& out , const std::unique_ptr<AbstractProposal> &rhs) 
{
  out << rhs->name <<  " primarily modifying " ; 
  for(auto &r : rhs->primVar)
    out << r.get() << ",\t"  ; 

  if(not rhs->secVar.empty() )
    {
      out << "\tand also modifying " ; 
      for(auto &r : rhs->secVar ) 
	out << r.get() << ",\t" ; 
    }
  return out; 
}


AbstractProposal::AbstractProposal( const AbstractProposal& rhs)
  : name(rhs.name)
  , sctr(rhs.sctr)
  , category(rhs.category)
  , relativeWeight(rhs.relativeWeight)
{
  for(auto &v : rhs.primVar)
    primVar.emplace_back(v->clone()); 
  for(auto &v : rhs.secVar)
    secVar.emplace_back(v->clone()); 
} 


