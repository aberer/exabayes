#include "AbstractProposal.hpp"


vector<RandomVariable*> AbstractProposal::getPrimVar() const
{
  vector<RandomVariable*> result; 
  for(auto &v : primVar)
    result.push_back(v.get()); 
  return result; 
}

 
vector<RandomVariable*> AbstractProposal::getSecVar() const 
{
  vector<RandomVariable*> result; 
  for(auto &v : secVar)
    result.push_back(v.get()); 
  return result; 
}


/**
   @brief valToAdd must not be on the log-scale 
 */ 
void AbstractProposal::updateHastings(double &hastings, double valToAdd, string whoDoneIt) 
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
