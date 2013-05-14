#include "BranchCollapser.hpp"
#include "eval.h"


BranchCollapser::BranchCollapser(double _relativeProbability)
{
  this->relativeProbability = _relativeProbability; //  the constructor argument should not have the exact same name as the member variable

  name = "BranchCollapser"; 	
  category = BRANCH_LENGTHS ; 	// check out categoryType.h
}


void BranchCollapser::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  branch b = rand.drawBranchUniform(traln); 

  modifiedBranch = b ; 
  nodeptr p = findNodeFromBranch(traln.getTr(),b) ; 
  modifiedBranch.length[0] = p->z[0];       
  
  auto brPr = prior.getBranchLengthPrior();       

  if(traln.isCollapsed(b))
    {      
      vector<double> zNews =  prior.drawFromPriorByCategory(BRANCH_LENGTHS, rand);
      assert(zNews.size() == 1); 
      b.length[0] = branchLengthToInternal(traln.getTr(), zNews[0]); 
      traln.setBranchLengthUnsafe(b); 
      updateHastings(hastings,  1 / exp(brPr->getLogProb(zNews)), "branchCollapser"); 

    }
  else
    {      
      nodeptr p = findNodeFromBranch(traln.getTr(), b);       
      double realZ = branchLengthToReal(traln.getTr(), p->z[0]); 
      vector<double> tmp = {realZ}; 
      updateHastings(hastings, exp(brPr->getLogProb(tmp)) , "branchCollapser") ;
      traln.collapseBranch(b); 

    }  

  prior.updateBranchLength(branchLengthToReal(traln.getTr(), modifiedBranch.length[0]) ,
			   branchLengthToReal(traln.getTr(), p->z[0] ) );
}


void BranchCollapser::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = findNodeFromBranch(traln.getTr(), modifiedBranch); 
  evaluateGenericWrapper(traln, p, FALSE);
}



void BranchCollapser::resetState(TreeAln &traln, PriorBelief &prior)
{  
  nodeptr p = findNodeFromBranch(traln.getTr(), modifiedBranch) ; 

  prior.updateBranchLength(branchLengthToReal(traln.getTr(), p->z[0] ) , 
			   branchLengthToReal(traln.getTr(), modifiedBranch.length[0]) );

  traln.setBranchLengthUnsafe(modifiedBranch);
}


AbstractProposal* BranchCollapser::clone() const 
{
  return new BranchCollapser(relativeProbability);
}

