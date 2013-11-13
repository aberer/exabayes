#include "TreeResource.hpp"
#include "TreeAln.hpp"
#include <cstring>


TreeResource::TreeResource(const TreeAln * tralnPtr)
  : _tralnPtr(tralnPtr)
{
}


std::vector<std::string> TreeResource::getTaxonNames(nat numTax) 
{
  auto taxa = _tralnPtr->getTaxa() ; 
  return taxa;
}

 
void TreeResource::fillAliasWgt(int *pos, nat length) 
{
}

   
std::tuple<int,int,double,int> TreeResource::getGlobalInfo()
{
  auto& tr = _tralnPtr->getTrHandle();
  auto numPart = _tralnPtr->getNumberOfPartitions(); 
  return std::make_tuple(
			 tr.mxtips,
			 numPart, 
			 tr.gapyness, 
			 tr.originalCrunchedLength
			 );
}


void TreeResource::fillPartition(pInfo &partition, nat model)  
{
  auto& rhsPart  = _tralnPtr->getPartition(model); 
  
  partition.states = rhsPart.states; 
  partition.maxTipStates = rhsPart.maxTipStates; 
  partition.lower = rhsPart.lower; 
  partition.upper = rhsPart.upper; 
  partition.width = rhsPart.width; 
  partition.dataType = rhsPart.dataType; 
  partition.protModels = rhsPart.protModels; 
  partition.protFreqs = rhsPart.protFreqs; 
  partition.nonGTR = rhsPart.nonGTR; 

  partition.partitionName = strdup(rhsPart.partitionName);
}


void TreeResource::fillAlnPart(unsigned char* ptr, nat length, nat &ctr) 
{
  assert(0);
  ++ctr; 
}

 
std::tuple<parsimonyNumber*, nat> TreeResource::fillParsVect(nat numNodes, nat states, nat model) 
{
  nat len = 0; 
  parsimonyNumber* ptr; 

  auto& partition  = _tralnPtr->getPartition(model);
  nat theseStates = partition.states; 
  assert(theseStates == states); 
  len = partition.parsimonyLength; 
  nat numBytes = numNodes * len * states; 
  // tout << "model=" << model << "\tcopying " << numBytes << std::endl; 
  ptr = (parsimonyNumber*)exa_malloc_aligned( numBytes * sizeof(parsimonyNumber));
  memcpy(ptr, partition.parsVect,numBytes * sizeof(parsimonyNumber));

  return std::make_tuple(ptr, len) ; 
}


void TreeResource::initWeightsAndAln(TreeAln &traln)  
{
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
    {
      auto &partitionLhs = traln.getPartition(i); 
      auto &partitionRhs = _tralnPtr->getPartition(i); 
      nat width = partitionRhs.width; 
      
      std::copy(partitionRhs.wgt, partitionRhs.wgt + width, partitionLhs.wgt); 
      
      for(nat j = 1; j < traln.getNumberOfTaxa() + 1 ; ++j)
	std::copy(partitionRhs.yVector[j] , partitionRhs.yVector[j] + width ,partitionLhs.yVector[j]);
    }
} 


std::vector<double> TreeResource::getPartitionContributions(nat num) 
{
  auto result = std::vector<double>(num, 0);
#if HAVE_PLL == 0
  auto start = _tralnPtr->getTrHandle().partitionContributions; 
  std::copy(start, start + num ,result.begin()); 
#else 
  for(nat i = 0; i< num ; ++i)
    {
      auto &partition = _tralnPtr->getPartition(i);
      assert(partition.partitionContribution > 0 ); 
      result.at(i) = partition.partitionContribution; 
    }
#endif
  return result; 
} 




