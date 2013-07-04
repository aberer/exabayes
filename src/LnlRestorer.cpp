#include <cstring>

#include <cassert>

#include "LnlRestorer.hpp" 


#include "Branch.hpp"
#include "TreeAln.hpp" 
#include "GlobalVariables.hpp"

#define GAMMA_CATS 4 


LnlRestorer::LnlRestorer(TreeAln& traln)
  : numPart(traln.getNumberOfPartitions())
  , numTax(traln.getTr()->mxtips)
  , wasSwitched( 2 * traln.getTr()->mxtips, false)
  , orientation( traln.getTr()->mxtips, 0)
  , partitionScaler(numPart, std::vector<nat> (2 * traln.getTr()->mxtips, 0) )
{
  tree *tr = traln.getTr(); 

  reserveArrays = (double***)exa_calloc(numPart, sizeof(double**)); 
  // partitionScaler = (nat**)exa_calloc(numPart, sizeof(nat*)); 
  // partitionScaler()
  
  partitionLnl = std::vector<double>( numPart);
  
  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition(i);
      reserveArrays[i] = (double**)exa_calloc(tr->mxtips, sizeof(double*)); 

#if HAVE_PLL != 0		/* TODO that's hacky */
      int length = partition->upper - partition->lower; 
#else 
      int length = partition->width; 
#endif

      for(int j = 0; j < numTax-2; ++j)
	{
	  reserveArrays[i][j] = (double*)exa_malloc_aligned(length *  partition->states * GAMMA_CATS  * sizeof(double)); // TODO not aligned? 
	  memset(reserveArrays[i][j], 0, length * partition->states * GAMMA_CATS * sizeof(double)); 
	}

      // partitionScaler[i] = (nat*)exa_calloc(2 * tr->mxtips , sizeof(nat)); 
    }  
}




LnlRestorer::~LnlRestorer()  
{
  for(int i = 0; i <  numPart; ++i)
    {
#if 0      
      // TODO ! 
      for(int j = 0; j < tr->mxtips; ++j)
	exa_free(reserveArrays[i][j]);
#endif

      exa_free(reserveArrays[i]); 
      
      // exa_free(partitionScaler[i]); 
    }
  exa_free(reserveArrays); 
  // exa_free(partitionScaler); 
  // exa_free(wasSwitched); 
}



/**
   @brief loads the saved orientation of the x-vectors
 */
void LnlRestorer::loadOrientation(TreeAln &traln)
{
  tree *tr = traln.getTr(); 
  int ctr = 0; 
  for(int i = tr->mxtips+1; i < 2 * tr->mxtips-1; ++i)
    {
      int val = orientation[ctr]; 
      Branch b = Branch(tr->nodep[i]->number, val ); 
      b.exists(traln); 
      
      nodeptr q = b.findNodePtr(traln);
      assert(q->back->number == val); 
      q->x = 1; q->next->x = 0; q->next->next->x = 0;
      
      ctr++; 
    } 
}


void LnlRestorer::storeOrientation(const TreeAln &traln)
{
  auto *tr = traln.getTr(); 
  int ctr = 0; 
  for(int i = tr->mxtips+1; i < 2 * tr->mxtips-1; ++i)
    {	
      nodeptr p = tr->nodep[i]; 
      if(p->x)
	orientation [ ctr] = p->back->number; 
      else if(p->next->x)
	orientation [ ctr] = p->next->back->number; 
      else if(p->next->next->x)
	orientation [ ctr] = p->next->next->back->number; 
      else 
	assert(0); 
      ctr++; 
    }
}



void LnlRestorer::swapArray(TreeAln &traln, int nodeNumber, int model)
{
  tree *tr = traln.getTr();   
  assert(not isTip(nodeNumber, tr->mxtips)); 
  
  if(model == ALL_MODELS)
    {
#ifdef DEBUG_ARRAY_SWAP      
      cout << "swapped array for node " << nodeNumber <<   " and all  models "; 
#endif

      for(int i = 0; i < numPart; ++i)
	{
	  pInfo *partition = traln.getPartition( i);
	  int posInArray = nodeNumber -( tr->mxtips + 1); 	  
	  double *&a =  reserveArrays[i][posInArray], 
	    *&b  = partition->xVector[posInArray] ; 
#ifdef DEBUG_ARRAY_SWAP	  
	  cout << a << "," << b << "\t"; 
#endif
	  if(not b )
	    return; 
	  std::swap(a,b); 	  
	}
#ifdef DEBUG_ARRAY_SWAP
      cout << endl; 
#endif
    }
  else 
    {
      pInfo *partition = traln.getPartition( model); 
      int posInArray = nodeNumber- (tr->mxtips+1); 
      double*& a =  reserveArrays[model][posInArray], 
	*&b  = partition->xVector[posInArray]; 

#ifdef  DEBUG_ARRAY_SWAP
      cout << "swapped array for node " << nodeNumber <<  " and model " << model << ": " << a   << "," << b  << endl; 
#endif
      if(not b )
	return; 
      std::swap(a,b );     
    }
}


void LnlRestorer::restoreArrays(TreeAln& traln)
{
#ifdef DEBUG_ARRAY_SWAP   
  cout << "RESTORE for model" << modelEvaluated << endl; 
#endif
  tree *tr = traln.getTr(); 
  // switch arrays 
  for(int i = 0; i < 2 * tr->mxtips; ++i)
    {
      if(wasSwitched[i])
	swapArray(traln,i , modelEvaluated); 
    }

  loadOrientation(traln);


  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition( i); 
      for(nat j = 0; j < 2 * traln.getNumberOfTaxa(); ++j)
	partition->globalScaler[j] =  partitionScaler[i][j]; 
    }

  tr->likelihood = prevLnl; 

  traln.setPartitionLnls(partitionLnl); 
}


void LnlRestorer::traverseAndSwitchIfNecessary(TreeAln &traln, nodeptr virtualRoot, int model, bool fullTraversal)
{  
  this->modelEvaluated = model; 
  tree *tr = traln.getTr(); 

  if(isTip(virtualRoot->number, tr->mxtips))
    return; 

  bool incorrect = not virtualRoot->x; 
  if( ( incorrect
	|| fullTraversal )      
      && not wasSwitched[virtualRoot->number])
    {
#ifdef DEBUG_ARRAY_SWAP
      cout << "incorr, unseen " << virtualRoot->number << endl; 
#endif
      wasSwitched[virtualRoot->number] = true; 
      
      swapArray(traln, virtualRoot->number, model); 
    }
  else if (incorrect)
    {
#ifdef DEBUG_ARRAY_SWAP
      
      cout << "incorr, seen " <<  virtualRoot->number  << endl; 
#endif
    }
#ifdef DEBUG_ARRAY_SWAP
  else
    cout << "corr "<<  virtualRoot->number << endl; 
#endif

  if(incorrect || fullTraversal)
    {
      traverseAndSwitchIfNecessary(traln, virtualRoot->next->back, model, fullTraversal); 
      traverseAndSwitchIfNecessary(traln, virtualRoot->next->next->back, model, fullTraversal); 
    }
}


void LnlRestorer::resetRestorer(const TreeAln &traln)
{
#ifdef DEBUG_ARRAY_SWAP
  cout << "RESETTING RESTORER" << endl; 
#endif

  auto *tr = traln.getTr();
  for(int i = 0 ; i < 2 * tr->mxtips; ++i)
    wasSwitched[i] = false;
  storeOrientation(traln); 
  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition( i); 
      for(nat j = 0; j < traln.getNumberOfTaxa(); ++j)
	partitionScaler[i][j] = partition->globalScaler[j]; 
    }
  modelEvaluated = ALL_MODELS; 
  prevLnl = tr->likelihood;   

  partitionLnl = traln.getPartitionLnls();
}

