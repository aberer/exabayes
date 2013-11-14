#include <cstring>
#include <algorithm>
#include <cassert>

#include "ArrayRestorer.hpp" 
#include "Branch.hpp"
#include "GlobalVariables.hpp"

#include "FlagType.hpp"

ArrayRestorer::ArrayRestorer(const TreeAln& traln, bool _cacheTipTip, bool _cacheTipInner)
  : restoresGapVector(false)
  , arrayOrientation(traln)
  , cacheTipTip(_cacheTipTip)
  , cacheTipInner(_cacheTipInner)
{
  bool useSEV = ( traln.getMode() & RunModes::MEMORY_SEV ) != RunModes::NOTHING; 

  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    partitionLikelihoods.emplace_back(traln,i, useSEV); 
}


void ArrayRestorer::restoreSomePartitions(TreeAln &traln, const std::vector<bool> &partitions, ArrayOrientation &evalOrientation)
{
  for(nat model = 0; model <  traln.getNumberOfPartitions() ;++model)
    {
      if(not partitions[model])
	continue; 
      nat partitionIndex = model; 
      // tout << "restoring partition "  << model << std::endl; 
      
      nat ctr = 0; 
      nat lastNode = traln.getNumberOfNodes() + 1; 
      for(nat i = traln.getNumberOfTaxa() + 1  ; i < lastNode ; ++i)
	{
	  if(partitionLikelihoods[partitionIndex].isCached.at(ctr))
	    {
	      uncache(traln, i, partitionIndex, evalOrientation); 
	    }
	  ++ctr; 
	}

      // restore the partition scaler 
      auto& partition = traln.getPartition( partitionIndex);       

      if(restoresGapVector)
	{
	  std::copy(partitionLikelihoods.at(partitionIndex).gapColumn.begin(), 
		    partitionLikelihoods.at(partitionIndex).gapColumn.end(), 
		    partition.gapColumn ); 
	}
    }
}


void ArrayRestorer::cache( TreeAln &traln, nat nodeNumber, nat partitionId, const ArrayOrientation &curOrient )
{
  // tout << "caching node=" << nodeNumber << ", partition=" << partitionId << std::endl; 

  auto id = nodeNumber - ( traln.getNumberOfTaxa() + 1 ) ; 

  if(partitionLikelihoods[partitionId].cachedArrays.at(id) != nullptr)
    {
      assert(0);
      exa_free(partitionLikelihoods[partitionId].cachedArrays.at(id)); 
    }

  auto arrayAndLength = removeArray(traln, nodeNumber, partitionId); 

  partitionLikelihoods[partitionId].cachedArrays[id] = arrayAndLength.first; 
  partitionLikelihoods[partitionId].lengths[id] = arrayAndLength.second; 
  
  auto& partition = traln.getPartition(partitionId); 
  partitionLikelihoods[partitionId].scaler[nodeNumber] = partition.globalScaler[nodeNumber]; 

  if(restoresGapVector)
    {
      auto& partition = traln.getPartition(partitionId); 
      auto vec = partition.gapVector + nodeNumber * partition.gapVectorLength; 
      auto iter = partitionLikelihoods[partitionId].gapVector.begin() + id * partition.gapVectorLength ; 
      std::copy(vec , vec + partition.gapVectorLength, iter ); 
    }
}


void ArrayRestorer::destroyAndForget(TreeAln &traln, nat nodeNumber, nat partitionId )
{
  auto id = nodeNumber - (traln.getNumberOfTaxa() + 1); 
  auto arrayAndLength = removeArray(traln, nodeNumber, partitionId);
  exa_free(arrayAndLength.first); 
  
  auto &elem  = partitionLikelihoods[partitionId];
  
  elem.cachedArrays[id] = nullptr; 
  elem.lengths[id] = 0; 

  arrayOrientation.setInvalid(partitionId, id);
}


void ArrayRestorer::uncache(TreeAln &traln, nat nodeNumber, nat partitionId, ArrayOrientation &curOrient )
{
  // tout << "uncaching node=" << nodeNumber << ", partition=" << partitionId << std::endl; 

  auto& partition = traln.getPartition(partitionId); 
  auto id = nodeNumber - ( traln.getNumberOfTaxa() + 1) ; 

  auto &backup = partitionLikelihoods[partitionId]; 
  auto arrayAndLength = std::make_pair(backup.cachedArrays.at(id) , 
				       backup.lengths.at(id)); 

  if(partition.xVector[id] != NULL)
    {
      exa_free(partition.xVector[id]); 
      partition.xVector[id] = NULL; 
      partition.xSpaceVector[id] = 0; 
    }

  insertArray(traln, nodeNumber, partitionId, arrayAndLength); 
  
  backup.cachedArrays.at(id) = nullptr; 
  backup.lengths.at(id) = 0; 
  
  // backup orientation 
  auto tmp = arrayOrientation.getOrientation(partitionId, id);
  curOrient.setOrientation(partitionId, id, tmp);

  partition.globalScaler[nodeNumber] = backup.scaler[nodeNumber]  ; 
  
  if(restoresGapVector)
    {
      auto& partition = traln.getPartition(partitionId); 
      auto vec = partition.gapVector + nodeNumber * partition.gapVectorLength; 
      auto iter = backup.gapVector.begin() + id * partition.gapVectorLength; 
      std::copy(iter, iter + partition.gapVectorLength, vec); 
    }
}


void ArrayRestorer::traverseAndCache(TreeAln &traln, nodeptr virtualRoot, nat model, ArrayOrientation& curOrient )
{  
  if(traln.isTipNode(virtualRoot))
    return; 

  nat index = virtualRoot->number - 1 - traln.getNumberOfTaxa(); 
  bool incorrect = not curOrient.isCorrect(model, index, virtualRoot->back->number);

  auto pnb = virtualRoot->next->back, 
    pnnb = virtualRoot->next->next->back;

  if( incorrect && not partitionLikelihoods[model].isCached[index] ) 
    {
      bool tiptip = traln.isTipNode(pnb) &&  traln.isTipNode(pnnb); 
      bool tipinner = ( not traln.isTipNode(pnb) )   !=  (   not traln.isTipNode(pnnb) ) ; 
      
      bool innerinner = not tiptip && not tipinner;   

      partitionLikelihoods.at(model).isCached[index] = true; 
	  
      if( ( tiptip && cacheTipTip ) 
	  || ( tipinner && cacheTipInner ) 
	  || innerinner )
	{
	  cache(traln, virtualRoot->number , model, curOrient); 
	}
      else 
	{
	  destroyAndForget(traln,virtualRoot->number, model);
	}
    }

  if(incorrect)
    {
      traverseAndCache(traln, pnb, model, curOrient); 
      traverseAndCache(traln, pnnb, model, curOrient); 
    }
}


void ArrayRestorer::resetRestorer(const TreeAln &traln, ArrayOrientation &curOrient)
{
  for(auto &p : partitionLikelihoods)
    p.isCached = std::vector<bool>(traln.getNumberOfInnerNodes(), false); 
  
  arrayOrientation = curOrient; 

  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      auto& partition = traln.getPartition( i); 

      if(restoresGapVector)
	{
	  // problematic? 
	  std::copy(partition.gapColumn ,
		    partition.gapColumn + traln.getNumberOfTaxa() * partition.states * 4, 
		    partitionLikelihoods.at(i).gapColumn.begin() ); 
	}
    }
}


void ArrayRestorer::clearMemory()
{
  for(auto &partitionLikelihood : partitionLikelihoods)
    {
      for(auto &elem : partitionLikelihood.cachedArrays )
	{
	  if(elem != nullptr)
	    exa_free(elem); 
	  elem = nullptr; 
	}
    }
}


std::pair<double*,nat> ArrayRestorer::removeArray(TreeAln &traln, nat num, nat pid)
{
  auto& partition = traln.getPartition(pid);

  nat id = num - traln.getNumberOfTaxa() - 1; 

  double *array = partition.xVector[id];
  double length = partition.xSpaceVector[id]; 
  
  partition.xVector[id] = NULL; 
  partition.xSpaceVector[id] = 0; 

  return std::make_pair(array, length); 
}


void ArrayRestorer::insertArray(TreeAln &traln, nat num, nat pid, std::pair<double*,nat> toInsert)
{
  auto& partition = traln.getPartition(pid);
  nat id = num - traln.getNumberOfTaxa() - 1 ; 

  assert(partition.xVector[id] == NULL); 
  
  partition.xVector[id] = toInsert.first; 
  partition.xSpaceVector[id] = toInsert.second; 
}
