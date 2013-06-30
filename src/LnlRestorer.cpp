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
{
  tree *tr = traln.getTr(); 

  reserveArrays = (double***)exa_calloc(numPart, sizeof(double**)); 
  partitionScaler = (nat**)exa_calloc(numPart, sizeof(nat*)); 
  
  partitionLnl = vector<double>( numPart);
  
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

      partitionScaler[i] = (nat*)exa_calloc(2 * tr->mxtips , sizeof(nat)); 
    }  
  
  orientation = (int*) exa_calloc(tr->mxtips, sizeof(int)); 
  wasSwitched = (bool*)exa_calloc( 2 * tr->mxtips , sizeof(bool));
}




LnlRestorer::~LnlRestorer()  
{
  assert(0); 
  for(int i = 0; i <  numPart; ++i)
    {
#if 0      
      // TODO ! 
      for(int j = 0; j < tr->mxtips; ++j)
	exa_free(reserveArrays[i][j]);
#endif

      exa_free(reserveArrays[i]); 
      
      exa_free(partitionScaler[i]); 
    }
  exa_free(reserveArrays); 
  exa_free(partitionScaler); 
  exa_free(wasSwitched); 
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
      int *val = orientation + ctr ; 
      nodeptr p = tr->nodep[i]; 
      if(p->x)
	*val = p->back->number; 
      else if(p->next->x)
	*val = p->next->back->number; 
      else if(p->next->next->x)
	*val = p->next->next->back->number; 
      else 
	assert(0); 
      ctr++; 
    }
}



void LnlRestorer::swapArray(TreeAln &traln, int nodeNumber, int model)
{
  tree *tr = traln.getTr();   
  assert(NOT isTip(nodeNumber, tr->mxtips)); 
  
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
	  if(NOT b )
	    return; 
	  swap(a,b); 	  
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
      if(NOT b )
	return; 
      swap(a,b );     
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
      memcpy(partition->globalScaler, partitionScaler[i], sizeof(nat) * 2 * traln.getTr()->mxtips);     
    }

  tr->likelihood = prevLnl; 
  for(int i = 0; i < numPart; ++i)
    traln.accessPartitionLH(i) = partitionLnl[i]; 
}


void LnlRestorer::traverseAndSwitchIfNecessary(TreeAln &traln, nodeptr virtualRoot, int model, bool fullTraversal)
{  
  this->modelEvaluated = model; 
  tree *tr = traln.getTr(); 

  if(isTip(virtualRoot->number, tr->mxtips))
    return; 

  bool incorrect = NOT virtualRoot->x; 
  if( ( incorrect
	|| fullTraversal )      
      && NOT wasSwitched[virtualRoot->number])
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
      memcpy(partitionScaler[i], partition->globalScaler, sizeof(nat) * 2 * traln.getTr()->mxtips);     
    }
  modelEvaluated = ALL_MODELS; 
  prevLnl = tr->likelihood;   
  
  for(int i = 0; i < numPart; ++i)
    partitionLnl[i] = traln.accessPartitionLH(i); 
}

