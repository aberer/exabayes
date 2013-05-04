#include "ExtendedSPR.hpp"
#include "output.h"
#include "topology-utils.h"
#include "TreeAln.hpp"


ExtendedSPR::ExtendedSPR(Chain *_chain, double _relativeWeight, double _stopProb, double _multiplier)
  : chain(_chain), stopProb(_stopProb), multiplier(_multiplier)    
{
  this->relativeProbability = _relativeWeight; 
  this->name = "eSPR"; 
  category = TOPOLOGY; 
  ptype = E_SPR; 
  modifiedPath = NULL;   
  createStack(&(this->modifiedPath));
}



ExtendedSPR::~ExtendedSPR()
{
  freeStack(&modifiedPath); 
}


/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
void ExtendedSPR::applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand)
{
  debug_printTree(chain);

  clearStack(modifiedPath); 
  drawPathForESPR(chain,modifiedPath,stopProb); 

  saveBranchLengthsPath(chain, modifiedPath); 

  applyPathAsESPR(chain->traln, modifiedPath);

#ifdef ESPR_MULTIPLY_BL
  multiplyAlongBranchESPR(chain, modifiedPath, multiplier);
#endif

  debug_checkTreeConsistency(chain->traln->getTr()); 
}



void ExtendedSPR::resetState(TreeAln &traln, PriorManager &prior )
{
  resetAlongPathForESPR (&traln, modifiedPath);   

  /* TODO resetAlong... should be able to correctly restore branch lengths    */
  restoreBranchLengthsPath(&traln, modifiedPath); 

  debug_checkTreeConsistency((&traln)->getTr());

  debug_printTree(chain); 
}




void ExtendedSPR::evaluateProposal(TreeAln &traln, PriorManager &prior)
{  
  tree *tr = traln.getTr(); 
  path *rPath = modifiedPath; 

  branch futureRoot = getThirdBranch(tr, rPath->content[0], rPath->content[1]); 
  
  /* evaluate at root of inserted subtree */
  nodeptr toEval = findNodeFromBranch(tr, futureRoot); /* dangerous */

  destroyOrientationAlongPath(tr, rPath, toEval); 
  destroyOrientationAlongPath(tr,rPath, toEval->back);

  /* printf("evaluating at branch %d,%d\n", toEval->number, toEval->back->number);  */
  evaluateGenericWrapper(chain, toEval, FALSE);
}

