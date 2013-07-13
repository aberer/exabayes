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
  , wasSwitched( traln.getNumberOfNodes(), false)
  , orientation( traln.getTr()->mxtips, 0)
  , partitionScaler(numPart )
{
  reservePerPartition.resize(numPart);   
  partitionLnl = std::vector<double>( numPart);
  
  for(int i = 0; i < numPart; ++i)
    {
      partitionScaler[i] = std::vector<nat> (traln.getNumberOfNodes(), 0);
      auto *partition = traln.getPartition(i);      
      auto &reserve = reservePerPartition[i]; 
      
#if HAVE_PLL != 0		/* TODO that's hacky */
      int length = partition->upper - partition->lower; 
#else 
      int length = partition->width; 
#endif

      for(int j = 0; j < numTax-2; ++j)
	{ 
	  size_t numBytes = length *  partition->states * GAMMA_CATS  * sizeof(double);
	  auto array = (double*)exa_malloc_aligned(numBytes); 

	  memset(array, 0, numBytes); 
	  reserve.push_back(array);
	}      
    }  
}



LnlRestorer::~LnlRestorer()
{
  for(auto &elem : reservePerPartition)
    for(auto &e : elem)
      exa_free(e); 
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
      assert(b.exists(traln)); 
      
      nodeptr q = b.findNodePtr(traln);
      assert(q->back->number == val); 

      q->x = 1; 
      q->next->x = 0; 
      q->next->next->x = 0;
      
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
	orientation[ctr] = p->back->number; 
      else if(p->next->x)
	orientation[ctr] = p->next->back->number; 
      else if(p->next->next->x)
	orientation[ctr] = p->next->next->back->number; 
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
      std::cout << "swapped array for node " << nodeNumber <<   " and all  models; pos= " << nodeNumber - (tr->mxtips + 1 ) << std::endl; 
#endif

      for(int i = 0; i < numPart; ++i)
	{
	  pInfo *partition = traln.getPartition( i);
	  int posInArray = nodeNumber -( tr->mxtips + 1); 
	  double *&a = reservePerPartition.at(i).at(posInArray), 
	    *&b  = partition->xVector[posInArray] ; 

	  std::swap(a,b); 	  
	}
    }
  else 
    {
      auto partition = traln.getPartition( model); 
      auto posInArray = nodeNumber - (tr->mxtips+1); 

#ifdef  DEBUG_ARRAY_SWAP
      std::cout << "swapped array for node " << nodeNumber <<  " and model " << model << std::endl; 
#endif

      double *&a = reservePerPartition.at(model).at(posInArray),
	*&b = partition->xVector[posInArray]; 

      std::swap(a,b );     
    }  
}


void LnlRestorer::restoreArrays(TreeAln& traln)
{
#ifdef DEBUG_ARRAY_SWAP   
  std::cout << "RESTORE for model" << modelEvaluated << std::endl; 
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
      for(nat j = 0; j < traln.getNumberOfNodes(); ++j)
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
      wasSwitched[virtualRoot->number] = true;       
      swapArray(traln, virtualRoot->number, model); 
    }

  if(incorrect || fullTraversal)
    {
      traverseAndSwitchIfNecessary(traln, virtualRoot->next->back, model, fullTraversal); 
      traverseAndSwitchIfNecessary(traln, virtualRoot->next->next->back, model, fullTraversal); 
    }
}


void LnlRestorer::resetRestorer(const TreeAln &traln)
{
#ifdef DEBUG_ARRAY_SWAP
  std::cout << "RESETTING RESTORER" << std::endl; 
#endif

  auto *tr = traln.getTr();
  wasSwitched = std::vector<bool>(traln.getNumberOfNodes(), false); 

  storeOrientation(traln); 
  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition( i); 
      for(nat j = 0; j < traln.getNumberOfNodes() ; ++j )
	partitionScaler[i][j] = partition->globalScaler[j]; 
    }
  modelEvaluated = ALL_MODELS; 
  prevLnl = tr->likelihood;   

  partitionLnl = traln.getPartitionLnls();
}

