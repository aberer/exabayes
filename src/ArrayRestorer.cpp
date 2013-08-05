#include <cstring>
#include <algorithm>
#include <cassert>

#include "ArrayRestorer.hpp" 
#include "Branch.hpp"
#include "TreeAln.hpp" 
#include "GlobalVariables.hpp"

#define GAMMA_CATS 4 


ArrayRestorer::ArrayRestorer(const TreeAln& traln)
  : numPart(traln.getNumberOfPartitions())
  , numTax(traln.getNumberOfTaxa())
  , orientation( traln.getNumberOfTaxa(),  0)
  , partitionScaler(numPart )
{
  reservePerPartition.resize(numPart);   

  for(nat i = 0; i < numPart; ++i)
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
  
  std::vector<bool> wasSwitched( traln.getNumberOfInnerNodes(), false); 
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    wasSwitchedByPartition.push_back(wasSwitched); 
}


ArrayRestorer::~ArrayRestorer()
{
  for(auto &elem : reservePerPartition)
    for(auto &e : elem)
      exa_free(e); 
}


void ArrayRestorer::loadOrientation(TreeAln &traln)
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


void ArrayRestorer::storeOrientation(const TreeAln &traln)
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



void ArrayRestorer::swapArray(TreeAln &traln, int nodeNumber, std::vector<nat> model)
{
  tree *tr = traln.getTr();   
  nat numTax = traln.getNumberOfTaxa();
  assert(not isTip(nodeNumber, tr->mxtips)); 

  for(auto &p : model )
    {
      auto partition = traln.getPartition(p); 
      auto posInArray = nodeNumber - (numTax+1); 

#ifdef DEBUG_ARRAY_SWAP
      if(debugPrint)
	std::cout << "swapped array for node " << nodeNumber <<  " and model " << p << std::endl; 
#endif

      double *&a = reservePerPartition.at(p).at(posInArray),
	*&b = partition->xVector[posInArray]; 

      std::swap(a,b );     
    }  
}


void ArrayRestorer::restoreSomePartitions(TreeAln &traln, std::vector<nat> partitions)
{
  for(auto &partitionIndex : partitions)
    {
      for(nat i = 0 ; i < traln.getNumberOfInnerNodes() ; ++i)
	{
	  nat nodeNumber = i + 1 + traln.getNumberOfTaxa(); 
	  if(wasSwitchedByPartition[partitionIndex][i])
	    swapArray(traln, nodeNumber , {partitionIndex}); 
	}
      
      // forget that we evaluated this partition 
      modelsEvaluated.erase(partitionIndex); 

      // restore the partition scaler 
      pInfo *partition = traln.getPartition( partitionIndex);       
      memcpy(partition->globalScaler, &(partitionScaler[partitionIndex][0]), sizeof(nat) * traln.getNumberOfTaxa() * 2 );
    }
}


void ArrayRestorer::restoreArrays(TreeAln& traln)
{
  std::vector<nat> listOfEvaluated; 
  for(auto &p: modelsEvaluated)
    listOfEvaluated.push_back(p);
  restoreSomePartitions(traln, listOfEvaluated);  
  loadOrientation(traln);
}


void ArrayRestorer::toplevelSwitch(TreeAln &traln, Branch virtualRoot, std::vector<nat> models, bool fullTraversal)
{  
  modelsEvaluated.insert(models.begin(),models.end() ); 

  traverseAndSwitchIfNecessary(traln, virtualRoot.findNodePtr(traln), models, fullTraversal); 
  traverseAndSwitchIfNecessary(traln, virtualRoot.getInverted().findNodePtr(traln), models, fullTraversal); 
}


void ArrayRestorer::traverseAndSwitchIfNecessary(TreeAln &traln, nodeptr virtualRoot, std::vector<nat> models, bool fullTraversal)
{  
  
  modelsEvaluated.insert(models.begin(), models.end()); 
  
  if(traln.isTipNode(virtualRoot))
    return; 

  bool incorrect = not virtualRoot->x; 
  nat index = virtualRoot->number - 1 - traln.getNumberOfTaxa(); 

  
  for(auto &partition : models)
    {
      if( ( incorrect || fullTraversal )      
	  && not wasSwitchedByPartition[partition][index ])
	{
	  wasSwitchedByPartition[partition][index] = true;       
	  swapArray(traln, virtualRoot->number, {partition}); 
	}
    }

  if(incorrect || fullTraversal)
    {
      traverseAndSwitchIfNecessary(traln, virtualRoot->next->back, models, fullTraversal); 
      traverseAndSwitchIfNecessary(traln, virtualRoot->next->next->back, models, fullTraversal); 
    }
}


void ArrayRestorer::resetRestorer(const TreeAln &traln)
{
  if(debugPrint)
    std::cout << "RESETTING RESTORER" << std::endl; 

  modelsEvaluated.clear(); 

  std::vector<bool> wasSwitched(traln.getNumberOfInnerNodes() , false); 
  for(auto &elem : wasSwitchedByPartition)
    elem = wasSwitched; 
  
  storeOrientation(traln); 

  for(nat i = 0; i < numPart; ++i)
    {
      pInfo *partition = traln.getPartition( i); 
      partitionScaler[i].assign(partition->globalScaler, 
				partition->globalScaler + traln.getNumberOfTaxa() * 2 ); 
    }
}
