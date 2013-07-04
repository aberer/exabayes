#include "LikelihoodEvaluator.hpp"
#include "GlobalVariables.hpp"
#include "LnlRestorer.hpp"


LikelihoodEvaluator::LikelihoodEvaluator(shared_ptr<LnlRestorer> _restorer)
  : restorer(_restorer)
{
}


void LikelihoodEvaluator::exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL != 0
  evaluateGeneric(traln.getTr(), traln.getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(traln.getTr(), start, fullTraversal); 
#endif  
}


  // BAD
void LikelihoodEvaluator::evaluateFullNoBackup(TreeAln& traln)
{
#ifdef DEBUG_EVAL  
  cout << "conducting full evaluation, no backup created" << endl; 
#endif
  
  exa_evaluateGeneric(traln,traln.getTr()->start,TRUE );   
#ifdef DEBUG_LNL_VERIFY
  expensiveVerify(traln);
#endif
} 



// must be partial 
void LikelihoodEvaluator::evalSubtree(TreeAln  &traln, const Branch &evalBranch)   
{ 
  bool masked = false; 
  nodeptr p = evalBranch.findNodePtr(traln); 
  int numberToExecute = 0; 
  int modelToEval = ALL_MODELS; 
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      if(traln.accessExecModel(i))
	{
	  numberToExecute++; 
	  modelToEval = i; 
	}
    }

  assert(numberToExecute == 1 || numberToExecute == traln.getNumberOfPartitions()); 
  if(numberToExecute > 1)
    modelToEval = ALL_MODELS; 

#ifdef DEBUG_EVAL  
  cout << "newViewGenericWrapper on " << p->number  << " and model " <<  modelToEval << endl; 
#endif

  if(p->x)
    {
      tree *tr = traln.getTr();
      assert(not isTip(p->number, tr->mxtips)); 
      p->x = 0; 
      p->next->x = 1; 
    }
  restorer->traverseAndSwitchIfNecessary(traln, p, modelToEval, false); 
  coreEvalSubTree(traln,p,masked); // NEEDED
}


double LikelihoodEvaluator::evaluate(TreeAln &traln, const Branch &evalBranch, bool fullTraversal )  
{
#ifdef DEBUG_EVAL  
  cout << "evaluateGeneric at " << start->number << "/" << start->back->number << " with " << (fullTraversal ? "TRUE" : "FALSE" )  << endl; 
#endif

  int model  = ALL_MODELS; 
  nodeptr start = evalBranch.findNodePtr(traln) ; 

  restorer->traverseAndSwitchIfNecessary(traln, start, model, fullTraversal);
  restorer->traverseAndSwitchIfNecessary(traln, start->back, model, fullTraversal);
  
  exa_evaluateGeneric(traln,start,fullTraversal);   
#ifdef DEBUG_LNL_VERIFY
  if(verifyLnl)
    expensiveVerify(traln);
#endif
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
  auto p = branch.findNodePtr(traln ); 

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
      restorer->traverseAndSwitchIfNecessary(traln,root.findNodePtr(traln), p, true ); 
      restorer->traverseAndSwitchIfNecessary(traln,root.getInverted().findNodePtr(traln), p, true); 
    }

  exa_evaluateGeneric(traln, root.findNodePtr(traln),  FALSE ); 
  
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

void LikelihoodEvaluator::coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked)
{
#if HAVE_PLL != 0
  newviewGeneric(traln.getTr(), traln.getPartitionsPtr(), p, masked); 
#else 
  newviewGeneric(traln.getTr(), p, masked); 
#endif 
}


#ifdef DEBUG_LNL_VERIFY
void LikelihoodEvaluator::expensiveVerify(TreeAln &traln)
{
  double toVerify = traln.getTr()->likelihood; 

  *debugTraln = traln; 

  // LikelihoodEvaluator eval; 
  Branch root ; 
  findVirtualRoot(traln, root); 

  nodeptr
    p = root.findNodePtr(*debugTraln ); 

  disorientTree(*debugTraln, root); 
  exa_evaluateGeneric(*debugTraln, p , FALSE); 

  double verifiedLnl =  debugTraln->getTr()->likelihood; 

  if(fabs (verifiedLnl - toVerify ) > ACCEPTED_LIKELIHOOD_EPS)
    {
      tout << "WARNING: found in expensive evaluation: likelihood difference is " 
	   << setprecision(8) <<   fabs (verifiedLnl - toVerify )
	   << " (with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")" << endl; 

      tout << "current tree: " << traln << endl; 
      tout << "help tree: " <<  *debugTraln << endl; 	  

      evaluateFullNoBackup(traln);
      tout << "full evaluation on original tree yields: "  << traln.getTr()->likelihood << endl;       
    }  
  assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS);   
}
#endif


#ifdef DEBUG_LNL_VERIFY
void LikelihoodEvaluator::setDebugTraln(shared_ptr<TreeAln> _debugTraln)
{
  verifyLnl = true; 
  debugTraln = _debugTraln; 
}
#endif





