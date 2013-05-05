#include "ExtendedTBR.hpp"
#include "branch.h"


ExtendedTBR::ExtendedTBR( Chain *_chain, double _relativeProbability, double _extensionProb, double _multiplier)
  : chain(_chain)
  ,  extensionProbability(_extensionProb)
  , multiplier(_multiplier)
{
  name = "eTBR"; 
  relativeProbability = _relativeProbability; 
  category = TOPOLOGY; 
  ptype = E_TBR; 
}



void ExtendedTBR::autotune()
{
  // cannot tune that 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorManager& prior, double &hastings, Randomness &rand)
{
  tree *tr = traln.getTr();

  branch bisectedBranch = {0,0}; 

  assert(tr->mxtips > 8 ); 

  nodeptr
    p1, p2;  

  do 
    {


      bisectedBranch = rand.drawInnerBranchUniform(traln); 
      p1 = findNodeFromBranch(tr, bisectedBranch); 
      p2 = findNodeFromBranch(tr, invertBranch( bisectedBranch)); 

      assert(NOT isTip(p1->number, tr->mxtips) && NOT isTip(p2->number, tr->mxtips)); 
      
      if(NOT ( isTip( p1->next->number, tr->mxtips) 
	       || isTip( p1->next->next->number, tr->mxtips)
	       || isTip( p2->next->number, tr->mxtips) 
	       || isTip( p2->next->next->number, tr->mxtips)))
	break; 
      
      
    } while( TRUE  ); 		// just because of the many nasty
				// variables otherwise

  // prune both nodes 
  double z1[NUM_BRANCHES],
    z2[NUM_BRANCHES]; 
  
  


}


void ExtendedTBR::evaluateProposal(TreeAln& traln, PriorManager& prior)
{

}


void ExtendedTBR::resetState(TreeAln &traln, PriorManager& prior)
{
  
}

