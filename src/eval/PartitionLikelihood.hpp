#ifndef _PARTITION_LIKELIHOOD 
#define _PARTITION_LIKELIHOOD 

#include <vector>
#include "model/TreeAln.hpp"

class PartitionLikelihood
{
public: 
  PartitionLikelihood(const TreeAln& traln, nat model, bool useSEV)
    : model(model)
    , cachedArrays(traln.getNumberOfInnerNodes()) 
    , lengths(traln.getNumberOfInnerNodes(), 0)
    , scaler( 2 * traln.getNumberOfTaxa() , 0)
    , isCached(traln.getNumberOfInnerNodes(),false)
  {
    if(useSEV)
      {
// #if 0
	auto& partition = traln.getPartition(model); 
	gapVector = std::vector<nat>(partition.getGapVectorLength() * 2 * traln.getNumberOfTaxa() , 0);
	gapColumn = std::vector<double>(traln.getNumberOfTaxa( ) * partition.getStates() * Partition::maxCategories ,0);
	// very MEH 

// #else 
// 	// TODO gap stuff 
// 	assert(0); 
// #endif
      }
  }


public: 
  nat model; 
  std::vector< double* > cachedArrays; 
  std::vector<size_t> lengths;
  std::vector<nat> scaler; 
  std::vector<bool> isCached; 
  std::vector<nat> gapVector;
  std::vector<double> gapColumn; 
}; 
#endif
