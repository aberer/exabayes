#include "TreeResource.hpp"
#include "TreeAln.hpp"
#include <cstring>


TreeResource::TreeResource(const TreeAln * tralnPtr)
  : _tralnPtr(tralnPtr)
{
}


std::vector<std::string> TreeResource::getTaxonNames(nat numTax) 
{
  return _tralnPtr->getTaxa();
}

 
void TreeResource::fillAliasWgt(int *pos, nat length) 
{
  auto &tr = _tralnPtr->getTrHandle();
  memcpy(pos, tr.aliaswgt, length * sizeof(int)); 
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
  
  // assert(0);
  
  ++ctr; 
}

 
void TreeResource::fillParsVect(parsimonyNumber*& ptr, size_t &len, nat mult, nat model)
{
  auto& partition  = _tralnPtr->getPartition(model);
  len = partition.parsimonyLength; 
  nat numBytes = mult * len; 
  ptr = (parsimonyNumber*)exa_malloc_aligned( numBytes * sizeof(parsimonyNumber));
  memcpy(ptr, partition.parsVect,numBytes * sizeof(parsimonyNumber));
}
