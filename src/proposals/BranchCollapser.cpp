#include "BranchCollapser.hpp"


BranchCollapser::BranchCollapser()
{
  name = "BranchCollapser"; 	
  category = Category::BRANCH_LENGTHS ; 	// check out categoryType.h

  relativeWeight = 0.0 ; 
}


void BranchCollapser::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  Branch b = traln.drawBranchUniform(rand); 

  modifiedBranch = b ; 
  nodeptr p = b.findNodePtr(traln) ; 
  modifiedBranch.setLength( p->z[0]); 

  auto brPr =   secVar[0]->getPrior();

  if(traln.isCollapsed(b))
    {    

      std::vector<double> zNews ; 
#if 0 
#if TODO
      zNews=  prior.drawFromPriorByCategory(BRANCH_LENGTHS, rand);
#endif
#endif

#if 0 
#if TODO 
      assert(zNews.size() == 1); 
      b.setLength( branchLengthToInternal(traln.getTr(), zNews[0])) ; 
      traln.setBranchLengthUnsafe(b); 
#endif
#endif
#if TODO 
      updateHastings(hastings,  1 / exp(brPr->getLogProb(zNews)), "branchCollapser"); 
#endif

    }
  else
    {      
      double realZ = b.getInterpretedLength(traln); 
      std::vector<double> tmp = {realZ}; 
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


void BranchCollapser::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, PriorBelief &prior) 
{
  evaluator.evaluate(traln,modifiedBranch, false); 
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

