#include "BranchCollapser.hpp"
#include "eval.h"





BranchCollapser::BranchCollapser()
{
  name = "BranchCollapser"; 	
  category = BRANCH_LENGTHS ; 	// check out categoryType.h

  relativeWeight = 0.0 ; 
}


void BranchCollapser::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  branch b = rand.drawBranchUniform(traln); 

  modifiedBranch = b ; 
  nodeptr p = findNodeFromBranch(traln.getTr(),b) ; 
  modifiedBranch.length[0] = p->z[0];       

  auto brPr =   secVar[0].getPrior();

  if(traln.isCollapsed(b))
    {    

      vector<double> zNews ; 
#if 0 
#if TODO
      zNews=  prior.drawFromPriorByCategory(BRANCH_LENGTHS, rand);
#endif
#endif

      assert(zNews.size() == 1); 
      b.length[0] = branchLengthToInternal(traln.getTr(), zNews[0]); 
      traln.setBranchLengthUnsafe(b); 
#if TODO 
      updateHastings(hastings,  1 / exp(brPr->getLogProb(zNews)), "branchCollapser"); 
#endif

    }
  else
    {      
      nodeptr p = findNodeFromBranch(traln.getTr(), b);       
      double realZ = branchLengthToReal(traln.getTr(), p->z[0]); 
      vector<double> tmp = {realZ}; 
#if TODO 
      updateHastings(hastings, exp(brPr->getLogProb(tmp)) , "branchCollapser") ;
#endif
      traln.collapseBranch(b); 

    }  

#if 0 
  prior.updateBranchLength(branchLengthToReal(traln.getTr(), modifiedBranch.length[0]) ,
			   branchLengthToReal(traln.getTr(), p->z[0] ) );
#endif
}


void BranchCollapser::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = findNodeFromBranch(traln.getTr(), modifiedBranch); 
  evaluateGenericWrapper(traln, p, FALSE);
}



void BranchCollapser::resetState(TreeAln &traln, PriorBelief &prior)
{  
#if 0 
#if TODO 
  nodeptr p = findNodeFromBranch(traln.getTr(), modifiedBranch) ; 
  prior.updateBranchLength(branchLengthToReal(traln.getTr(), p->z[0] ) , 
			   branchLengthToReal(traln.getTr(), modifiedBranch.length[0]) );
#endif
#endif

  traln.setBranchLengthUnsafe(modifiedBranch);
}


AbstractProposal* BranchCollapser::clone() const 
{
  return new BranchCollapser(*this);
}
