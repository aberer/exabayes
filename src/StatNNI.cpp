#include "StatNNI.hpp"
#include "eval.h"
#include "Path.hpp"
#include "output.h"



StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
{
  // this->relativeProbability = _weight; 
  this->name = "stNNI" ; 
  this->category = TOPOLOGY; 
  // ptype = ST_NNI; 
  
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
    if(not prior.believingInFixedBranchLengths() )
      {

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

#ifdef PRINT_MULT
	printf("%f*%f=%f\t%f*%f=%f\t%f*%f=%f\n", old1, m1, new1, old2, m2, new2, old3,m3,new3) ;
#endif

	prior.updateBranchLength(branchLengthToReal(traln.getTr(), old1), branchLengthToReal(traln.getTr(), new1));
	prior.updateBranchLength(branchLengthToReal(traln.getTr(), old2), branchLengthToReal(traln.getTr(), new2)); 
	prior.updateBranchLength(branchLengthToReal(traln.getTr(), old3), branchLengthToReal(traln.getTr(), new3)); 

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
    
    if(not prior.believingInFixedBranchLengths())
      {
	prior.updateBranchLength(branchLengthToReal(traln.getTr(), between->z[0]), branchLengthToReal(traln.getTr(), chosenBranch.length[0]));
	prior.updateBranchLength(branchLengthToReal(traln.getTr(), r->z[0]), branchLengthToReal(traln.getTr(), b.length[0]));
	prior.updateBranchLength(branchLengthToReal(traln.getTr(), q->z[0]), branchLengthToReal(traln.getTr(), a.length[0]));
      }

    traln.clipNode(between, between->back, chosenBranch.length[0]); 
    traln.clipNode(r, qBack, b.length[0]);  /* r->z */
    traln.clipNode(q, rBack, a.length[0]) ; /* q->z */
  }

  debug_checkTreeConsistency(traln.getTr());

  path.clear();
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( multiplier);
}
