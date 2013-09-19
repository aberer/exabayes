#include "LikelihoodEvaluator.hpp"
#include "GlobalVariables.hpp"
#include "ArrayRestorer.hpp"
#include "Branch.hpp"


void LikelihoodEvaluator::exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
  if(fullTraversal)
    tout << "FULL TRAVERSAL" << std::endl; 

#if HAVE_PLL != 0
  evaluateGeneric(traln.getTr(), traln.getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(traln.getTr(), start, fullTraversal); 
#endif  
}

double LikelihoodEvaluator::evaluatePartitions(TreeAln &traln, const std::vector<nat>& partitions, bool fullTraversal)
{
  auto root = findVirtualRoot(traln);
  // notice: cannot overload a virtual function?  
  return evaluatePartitionsWithRoot(traln, root, partitions,fullTraversal); 
}


void LikelihoodEvaluator::disorientTree(TreeAln &traln, const BranchPlain& root) 
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


void LikelihoodEvaluator::disorientSubtree(TreeAln &traln, const BranchPlain &branch) 
{  
  auto p = branch.findNodePtr(traln ); 

  disorientNode(p); 

  if(not traln.isTipNode(p))
    {
      disorientSubtree(traln, BranchPlain(p->next->back->number, p->number)); 
      disorientSubtree(traln, BranchPlain(p->next->next->back->number, p->number)); 
    }
}


// currently a bit expensive  
BranchPlain LikelihoodEvaluator::findVirtualRoot(const TreeAln &traln) const 
{  
  auto tr = traln.getTr(); 

  BranchPlain root; 
  for(int i = tr->mxtips +1 ; i < 2* tr->mxtips-1 ; ++i)
    {
      nodeptr
	p = tr->nodep[i],
	q = p;       
      do 
	{
	  auto newRoot = BranchPlain(q->number, q->back->number); 
	  if(q->x && q->back->x && not root.equalsUndirected(newRoot))
	    {
	      if(root.getPrimNode() != 0 )
		{
#if 0 
		  std::cout << "root already taken! " << root << " now at " << newRoot << std::endl; 
#else 
		  assert(0); 
#endif
		}
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
#if 0 
	      std::cout << "previous root was " << root << " now at " << Branch(p->number , p->back->number)  << std::endl; 	      
#else 
	      assert(0); 
#endif
	    }
	  assert(root.getPrimNode() == 0); 
	  root = BranchPlain(p->number, p->back->number); 
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

  auto root = findVirtualRoot(traln); 

  auto
    p = root.findNodePtr(*debugTraln ); 

  disorientTree(*debugTraln, root); 
  exa_evaluateGeneric(*debugTraln, p , FALSE); 

  double verifiedLnl =  debugTraln->getTr()->likelihood; 

  if(fabs (verifiedLnl - toVerify ) > ACCEPTED_LIKELIHOOD_EPS)
    {
      std::cerr  << "WARNING: found in expensive evaluation: likelihood difference is " 
	   <<  std::setprecision(8) <<   fabs (verifiedLnl - toVerify )
	   << " (with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")" << std::endl; 
      
      // what to print? 
      assert(0);
      // tout << "current tree: " << traln << std::endl; 
      // tout << "help tree: " <<  *debugTraln << std::endl; 	        
    }  
  assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS); 
}


#ifdef DEBUG_LNL_VERIFY
void LikelihoodEvaluator::setDebugTraln(std::shared_ptr<TreeAln> _debugTraln)
{
  verifyLnl = true; 
  debugTraln = _debugTraln; 
}
#endif


bool LikelihoodEvaluator::isDirty(nat partition,nat nodeId)  const 
{
  return dirty[partition][nodeId]; 
}


void LikelihoodEvaluator::setDirty(nat partition, nat nodeId, bool isDirty)
{
  dirty[partition][nodeId] = isDirty ; 
}
