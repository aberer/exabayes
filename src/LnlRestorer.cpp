#include <cstring>
#include <algorithm>
#include <cassert>

#include "LnlRestorer.hpp" 
#include "Branch.hpp"
#include "TreeAln.hpp" 
#include "GlobalVariables.hpp"

#define GAMMA_CATS 4 


LnlRestorer::LnlRestorer(const TreeAln& traln)
  : numPart(traln.getNumberOfPartitions())
  , numTax(traln.getNumberOfTaxa())
  , wasSwitched(traln.getNumberOfInnerNodes(), false)
  , orientation( traln.getNumberOfTaxa(),  0)
  , partitionScaler(numPart )
{
  reservePerPartition.resize(numPart);   

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

      for(nat j = 0; j < traln.getNumberOfInnerNodes(); ++j)
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
  nat numTax = traln.getNumberOfTaxa(); 

  for(nat i = 0; i < traln.getNumberOfInnerNodes(); ++i)
    {
      nat index = i + numTax + 1 ;  
      int val = orientation.at(i); 
      Branch b(tr->nodep[index]->number, val ); 
      assert(b.exists(traln)); 
      
      nodeptr q = b.findNodePtr(traln);
      assert(q->back->number == val); 
      
      q->x = 1; 
      q->next->x = 0; 
      q->next->next->x = 0;
    } 
}


void LnlRestorer::storeOrientation(const TreeAln &traln)
{
  auto *tr = traln.getTr(); 
  nat numTax = traln.getNumberOfTaxa(); 
  for(nat i = 0; i < traln.getNumberOfInnerNodes(); ++i)
    {	
      nat index = i + numTax + 1 ;
      nodeptr p = tr->nodep[index]; 
      if(p->x)
	orientation.at(i) = (nat)p->back->number; 
      else if(p->next->x)
	orientation.at(i) = (nat)p->next->back->number; 
      else if(p->next->next->x)
	orientation.at(i) = (nat)p->next->next->back->number; 
      else 
	assert(0); 

      assert(orientation.at(i) < traln.getNumberOfNodes() + 1  ); 
    }
}



void LnlRestorer::swapArray(TreeAln &traln, int nodeNumber, int model)
{
  tree *tr = traln.getTr();   
  nat numTax = traln.getNumberOfTaxa();
  assert(not isTip(nodeNumber, tr->mxtips)); 

  if(model == ALL_MODELS)
    {
#ifdef DEBUG_ARRAY_SWAP
      if(debugPrint)
	std::cout << "swapped array for node " << nodeNumber <<   " and all  models; pos= " << nodeNumber - (numTax + 1 ) << std::endl;       
#endif

      for(int i = 0; i < numPart; ++i)
	{
	  pInfo *partition = traln.getPartition( i);
	  int posInArray = nodeNumber -( numTax + 1); 
	  double *&a = reservePerPartition.at(i).at(posInArray), 
	    *&b  = partition->xVector[posInArray] ; 

	  std::swap(a,b); 	  
	}
    }
  else 
    {
      auto partition = traln.getPartition( model); 
      auto posInArray = nodeNumber - (numTax+1); 

#ifdef DEBUG_ARRAY_SWAP
      if(debugPrint)
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
  if(debugPrint)
    std::cout << "RESTORE for model" << modelEvaluated << std::endl; 
#endif

  // switch arrays 
  for(nat i = 0 ; i < traln.getNumberOfInnerNodes() ; ++i)
    {
      nat nodeNumber = i + 1 + traln.getNumberOfTaxa(); 
      if(wasSwitched[i])
	swapArray(traln, nodeNumber , modelEvaluated); 
    }

  loadOrientation(traln);


  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition( i);       
      memcpy(partition->globalScaler, &(partitionScaler[i][0]) , 
	     sizeof(nat) * traln.getNumberOfTaxa() * 2 );
    }

}


void LnlRestorer::traverseAndSwitchIfNecessary(TreeAln &traln, nodeptr virtualRoot, int model, bool fullTraversal)
{  
  this->modelEvaluated = model; 

  if(traln.isTipNode(virtualRoot))
    return; 

  bool incorrect = not virtualRoot->x; 
  nat index = virtualRoot->number - 1 - traln.getNumberOfTaxa(); 

  if( ( incorrect || fullTraversal )      
      && not wasSwitched[index ])
    {
      wasSwitched.at(index) = true;       
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
  if(debugPrint)
    std::cout << "RESETTING RESTORER" << std::endl; 

  wasSwitched = std::vector<bool>(traln.getNumberOfInnerNodes() , false); 
  storeOrientation(traln); 

  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition( i); 
      partitionScaler[i].assign(partition->globalScaler, 
				partition->globalScaler + traln.getNumberOfTaxa() * 2 ); 
    }

  modelEvaluated = ALL_MODELS; 
}


