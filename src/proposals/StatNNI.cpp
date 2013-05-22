#include "StatNNI.hpp"
#include "eval.h"
#include "Path.hpp"
#include "output.h"


double StatNNI::relativeWeight = 5 ;


StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "stNNI" ; 
  this->category = TOPOLOGY; 
}


void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{  
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 
  branch b = rand.drawInnerBranchUniform(traln); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  branch toBePruned = constructBranch(p->next->back->number, p->number ); 
  
  branch switchingBranch; 
  double r = rand.drawRandDouble01();
  double bl = 0; 
  if(r < 0.5)
    {
      switchingBranch = constructBranch(p->back->next->back->number, p->back->number); 
      bl = traln.getBranchLength( p->back->next->back,0); 
    }
  else 
    {
      switchingBranch = constructBranch(p->back->next->next->back->number, p->back->number); 
      bl = traln.getBranchLength( p->back->next->next->back,0); 
    }

#ifdef DEBUG_INFO
  cout << "stNNI: switching " <<  toBePruned <<  " with "  << switchingBranch << "\t common branch " << b << endl; 
#endif

  toBePruned.length[0] = traln.getBranchLength( p->next->back,0); 
  b.length[0] = traln.getBranchLength( p,0); 
  switchingBranch.length[0] = bl; 
  
  // Path *rememStack = pf->remembrance.modifiedPath; 
  path.clear(); 
  path.append( toBePruned);
  path.append( b); 
  path.append( switchingBranch);

  /* switch the branches  */
  {
    nodeptr q = findNodeFromBranch(tr, toBePruned),
      qBack = q->back; 
    
    nodeptr r = findNodeFromBranch(tr, switchingBranch),
      rBack = r->back; 
    
    
#ifdef NNI_MULTIPLY_BL
    assert(traln.getNumBranches() == 1 ); 
    if(traln.getBranchLengthsFixed()[0])
      {

    // check for fixed  

    // assert(NOT_IMPLEMENTED); 
    // if(not prior.believingInFixedBranchLengths() )
      // {

	/* DIRTY */
	double m1 = rand.drawMultiplier( multiplier), 
	  m2 =  rand.drawMultiplier(multiplier),
	  m3 =  rand.drawMultiplier(multiplier);
      

	double old1 = traln.getBranchLength( p,0),
	  old2 = traln.getBranchLength( r,0), 
	  old3 = traln.getBranchLength( q,0),
	  new1 = pow( old1,m1),
	  new2 = pow( old2 ,m2),
	  new3 = pow( old3 ,m3); 

	traln.setBranchLengthBounded( new1 ,0,p); 
	traln.setBranchLengthBounded( new2 ,0,r); 
	traln.setBranchLengthBounded( new3 ,0,q); 

	double realNew1  =traln.getBranchLength(p,0),
	  realNew2 = traln.getBranchLength(r,0),
	  realNew3 = traln.getBranchLength(q,0); 

	double realM1 = log(realNew1)/ log(old1),
	  realM2 = log(realNew2) / log(old2),
	  realM3 = log(realNew3) / log(old3); 

	updateHastings(hastings, realM1 * realM2 * realM3 , name); 
	
#ifdef PRINT_MULT
	printf("%f*%f=%f\t%f*%f=%f\t%f*%f=%f\n", old1, m1, new1, old2, m2, new2, old3,m3,new3) ;
#endif

	auto brPr = prior.getBranchLengthPrior(); 
	prior.updateBranchLengthPrior(traln, old1, realNew1, brPr);
	prior.updateBranchLengthPrior(traln, old2, realNew2, brPr);
	prior.updateBranchLengthPrior(traln, old3, realNew3, brPr);

      }
#endif

    traln.clipNode(p,p->back, p->z[0]); 
    traln.clipNode(r, qBack, r->z[0]);
    traln.clipNode(q, rBack, q->z[0]);    

  }
  
  debug_checkTreeConsistency(traln.getTr());
  /* TODO maybe multiply as well */
  
}




void StatNNI::evaluateProposal(TreeAln &traln, PriorBelief &prior)
{
  tree *tr = traln.getTr(); 

  branch b1 = path.at(0), // tobepruned 
    b2  = path.at(2);	// switchingbranch 
  
  branch exchangeBranch = constructBranch(b1.thatNode, b2.thatNode);
  
  nodeptr p = findNodeFromBranch(tr, exchangeBranch),
    q = findNodeFromBranch(tr, invertBranch(exchangeBranch)); 

  newViewGenericWrapper(traln, p, FALSE); 
  newViewGenericWrapper(traln, q, FALSE); 

  evaluateGenericWrapper(traln,p,FALSE ); 
}



void StatNNI::resetState(TreeAln &traln, PriorBelief &prior)  
{
  tree
    *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 

  branch b = path.at(0),
    chosenBranch =path.at(1),
    a = path.at(2); 

  swap(a.thatNode, b.thatNode); 


  /* reset the branches */
  {
    nodeptr q = findNodeFromBranch(tr, a),
      qBack = q->back; 
    nodeptr r = findNodeFromBranch(tr, b), 
      rBack = r->back; 
    
    nodeptr between = findNodeFromBranch(tr, chosenBranch ); 
    
    assert(numBranches == 1); 

    traln.clipNode(between, between->back, chosenBranch.length[0]); 
    traln.clipNode(r, qBack, b.length[0]);  /* r->z */
    traln.clipNode(q, rBack, a.length[0]) ; /* q->z */
  }

  debug_checkTreeConsistency(traln.getTr());

  path.clear();
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
