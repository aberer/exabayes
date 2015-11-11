#include "LikelihoodEvaluator.hpp"
#include "GlobalVariables.hpp"
#include "LnlRestorer.hpp"


void LikelihoodEvaluator::exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL != 0
  evaluateGeneric(traln.getTr(), traln.getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(traln.getTr(), start, fullTraversal); 
#endif  
}


void LikelihoodEvaluator::disorientTree(TreeAln &traln, const Branch &root) 
{
  disorientSubtree(traln,root);   
  disorientSubtree(traln, root.getInverted()); 
}


bool LikelihoodEvaluator::disorientNode( nodeptr p)
{
  bool result =  ( p->x == 1 ) ; 
  if(p->x)
    {
      p->x = 0; 
      p->next->x = 1; 
    }
  return result; 
}


void LikelihoodEvaluator::disorientSubtree(TreeAln &traln, const Branch &branch) 
{  
  auto p = branch.findNodePtr(traln ); 

  disorientNode(p); 

  if(not traln.isTipNode(p))
    {
      disorientSubtree(traln, Branch(p->next->back->number, p->number,0)); 
      disorientSubtree(traln, Branch(p->next->next->back->number, p->number,0)); 
    }
}


// currently a bit expensive  
Branch LikelihoodEvaluator::findVirtualRoot(const TreeAln &traln) const 
{  
  auto tr = traln.getTr(); 

  Branch root; 
  for(int i = tr->mxtips +1 ; i < 2* tr->mxtips-1 ; ++i)
    {
      nodeptr
	p = tr->nodep[i],
	q = p;       
      do 
	{
	  Branch newRoot(q->number, q->back->number); 
	  if(q->x && q->back->x && not root.equalsUndirected(newRoot))
	    {
	      if(root.getPrimNode() != 0 )
		std::cout << "root already taken! " << root << " now at " << newRoot << std::endl; 
	      assert(root.getPrimNode( )== 0 ) ;
	      root = newRoot; 
	    }
	  q = q->next; 
	} while(p != q); 
    }


  for(int i = 1; i < tr->mxtips+1; ++i)
    {
      nodeptr
	p = tr->nodep[i]; 
      if(p->back->x)
	{
	  if(root.getPrimNode() != 0)
	    {	      
	      std::cout << "previous root was " << root << " now at " << Branch(p->number , p->back->number)  << std::endl; 	      
	    }
	  assert(root.getPrimNode() == 0); 
	  root = Branch(p->number, p->back->number); 
	}
    }

  assert(root.getPrimNode( )!= 0); 

  return root; 
}


void LikelihoodEvaluator::coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked)
{
#if HAVE_PLL != 0
  newviewGeneric(traln.getTr(), traln.getPartitionsPtr(), p, masked); 
#else 
  newviewGeneric(traln.getTr(), p, masked); 
#endif 
}


void LikelihoodEvaluator::expensiveVerify(TreeAln &traln)
{
  double toVerify = traln.getTr()->likelihood; 

  debugTraln->copyModel(traln); 

  Branch root = findVirtualRoot(traln); 

  nodeptr
    p = root.findNodePtr(*debugTraln ); 

  disorientTree(*debugTraln, root); 
  exa_evaluateGeneric(*debugTraln, p , FALSE); 

  double verifiedLnl =  debugTraln->getTr()->likelihood; 

  if(fabs (verifiedLnl - toVerify ) > ACCEPTED_LIKELIHOOD_EPS)
    {
      tout << "WARNING: found in expensive evaluation: likelihood difference is " 
	   <<  std::setprecision(8) <<   fabs (verifiedLnl - toVerify )
	   << " (with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")" << std::endl; 

      tout << "current tree: " << traln << std::endl; 
      tout << "help tree: " <<  *debugTraln << std::endl; 	        
    }  
  assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS); 

  // tout << "VERIFIED " << traln.getTr()->likelihood << std::endl; 

  // TODO remove 
  // disorientTree(traln, root); 
  // exa_evaluateGeneric(traln, root.findNodePtr(traln) , FALSE);   
}


#ifdef DEBUG_LNL_VERIFY
void LikelihoodEvaluator::setDebugTraln(std::shared_ptr<TreeAln> _debugTraln)
{
  verifyLnl = true; 
  debugTraln = _debugTraln; 
}
#endif

