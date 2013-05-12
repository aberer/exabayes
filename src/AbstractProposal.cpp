#include "AbstractProposal.hpp"

void updateHastings(double &hastings, double valToAdd, string whoDoneIt) 
{
#ifdef DEBUG_HASTINGS  
  tout << whoDoneIt << "updates hastings " << hastings << " with " << valToAdd ; 
#endif

  hastings += log(valToAdd); 	// we are logarithmic now   

#ifdef DEBUG_HASTINGS
  tout << " => " << hastings << endl; 
#endif
}
