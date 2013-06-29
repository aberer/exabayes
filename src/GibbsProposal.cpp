#include "GibbsProposal.hpp"

#include "AbstractProposal.hpp"

const double GibbsProposal::EST_EXP = -0.87246; 
const double GibbsProposal::EST_FAC = 0.008067; 


/** 
    @brief draws a branch from the estimated posterior probability   
    
    internal length used as a starting point for newton-raphson
    optimization

    important: uses the current length for the hastings

    @param Branch branch -- contains the initial branch length
    
 */ 
double GibbsProposal::drawFromEsitmatedPosterior(Branch &branch, TreeAln &traln, Randomness& rand, double initVal,  int maxIter, double &hastings) 
{
  auto p = branch.findNodePtr(traln); 
  double initLength = branch.getLength(); 
  
  double nrD2 = 0; 
  branch.optimise(traln, nrD2, maxIter);

  double nrOpt = branch.getLength(); 
  double c = EST_FAC * pow(-nrD2, EST_EXP); 
  double beta = nrOpt + sqrt(pow(nrOpt,2) + 4 * c) / (2 * c ); 
  double alpha = nrOpt  * beta + 1 ; 
  
  double proposal = rand.drawRandGamma(alpha, beta);
  double prop = branch.getInternalLength(traln, proposal );   
  AbstractProposal::updateHastings(hastings, log(initLength) / log(prop), "theGibbs");
  branch.setLength(prop); 
  
  return 0 ; 
}
