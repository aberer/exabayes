#include "LikelihoodEvaluator.hpp"
#include "eval.h"
#include "GlobalVariables.hpp"
#include "LnlRestorer.hpp"

// TODO const correctness  
// TODO use move lnl restorer into this here 

// must be partial 
void LikelihoodEvaluator::evalSubtree(TreeAln  &traln, const Branch &evalBranch)   
{ 
  newViewGenericWrapper(traln, evalBranch.findNodeFromBranch(traln), false); 
}


double LikelihoodEvaluator::evaluate(TreeAln &traln, const Branch &evalBranch, bool fullTraversal )  
{
  evaluateGenericWrapper(traln, evalBranch.findNodeFromBranch(traln), fullTraversal ? TRUE : FALSE );  
  return traln.getTr()->likelihood;  
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
  auto p = branch.findNodeFromBranch(traln ); 

  disorientNode(p); 

  if(not traln.isTipNode(p))
    {
      disorientSubtree(traln, Branch(p->next->back->number, p->number,0)); 
      disorientSubtree(traln, Branch(p->next->next->back->number, p->number,0)); 
    }
}


// currently a bit expensive  
void LikelihoodEvaluator::findVirtualRoot(const TreeAln &traln, Branch &result) const 
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
		cout << "root already taken! " << root << " now at " << newRoot << endl; 
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
	      cout << "previous root was " << root << " now at " << Branch(p->number , p->back->number)  << endl; 	      
	    }
	  assert(root.getPrimNode() == 0); 
	  root = Branch(p->number, p->back->number); 
	}
    }

  assert(root.getPrimNode( )!= 0); 

  result = root; 
}


double LikelihoodEvaluator::evaluatePartitions( TreeAln &traln, const vector<nat>& partitions)    
{
  Branch root; 
  findVirtualRoot(traln,root); 
  
  auto tr = traln.getTr(); 
  nat numPart = traln.getNumberOfPartitions(); 

  vector<double> perPartitionLH; 
  for(nat i = 0; i < numPart; ++i )
    perPartitionLH.push_back(traln.accessPartitionLH(i));
  
  for(nat i = 0; i < numPart; ++i)
    traln.accessExecModel(i) = FALSE; 
  
  for(auto m : partitions)
    traln.accessExecModel(m) = TRUE; 

  disorientTree(traln, root); 
  
  for(auto p : partitions)
    {
      traln.getRestorer()->traverseAndSwitchIfNecessary(traln,root.findNodeFromBranch(traln), p, true ); 
      traln.getRestorer()->traverseAndSwitchIfNecessary(traln,root.getInverted().findNodeFromBranch(traln), p, true); 
    }

  exa_evaluateGeneric(traln, root.findNodeFromBranch(traln),  FALSE ); 
  
  for(auto m : partitions )
    perPartitionLH[m] = traln.accessPartitionLH(m); 

  for(nat i = 0; i < numPart; ++i )
    traln.accessPartitionLH(i) = perPartitionLH[i]; 

  tr->likelihood = 0; 
  for(nat i = 0; i < numPart; ++i)
    {
      tr->likelihood += traln.accessPartitionLH(i); 
      traln.accessExecModel(i) = TRUE; 
    }  
  
#ifdef DEBUG_LNL_VERIFY
  expensiveVerify(traln);   
#endif
  
  return tr->likelihood; 
}

