#include "ConsensusTree.hpp"
#include <map> 
#include <cassert> 
#include <algorithm>


ConsensusTree::ConsensusTree(std::vector<std::string> files, nat burnin  )
  : bipEx(files, false)
{
  bipEx.extractBips<false>(burnin); 
}


std::vector<Bipartition> ConsensusTree::getRefinedConsensusTree(const std::vector<Bipartition> &consensusBips, const std::vector<Bipartition> &minorityBips) const 
{
  auto result = consensusBips; 
  nat maxBipsNeeded = bipEx.getTaxa().size()  - 3 ; 
  nat maxElem = bipEx.getTaxa().size(); 
  
  for(auto minBip : minorityBips)
    {
      if(result.size() >= maxBipsNeeded)
	break; 

      auto isCompat = true; 
      for(auto &bip : result)
	{
	  isCompat &= bip.isCompatible(minBip, maxElem);
	  if(not isCompat)
	    break; 
	}
      
      if(isCompat)
	result.push_back(minBip); 
    }

  return result; 
}

					
std::string ConsensusTree::getConsensusTreeString(double threshold, bool isMRE)
{
  auto result = std::string{}; 
  const auto& bipHashes = bipEx.getBipartitionHashes();

  auto bip2Occ = std::unordered_map<Bipartition,nat>{}; 
  nat maxBip = bipEx.getTaxa().size() - 3; 

  nat totalTreesAdded = 0; 
  for(auto &bipHash : bipHashes)
    totalTreesAdded += bipHash.getTreesAdded(); 

  nat absThreshold = 0; 
  if(isMRE)
    absThreshold = nat(totalTreesAdded * .5);
  else 
    absThreshold = nat(totalTreesAdded * threshold); 

  assert(absThreshold > 0. ); 


  // determine unique bipartitions  
  for(auto &bipHash : bipHashes)
    {
      for(auto &bipOccPair : bipHash)
	{
	  auto &bip = std::get<0>(bipOccPair); 
	  if( bip2Occ.find(bip) == bip2Occ.end())
	    bip2Occ[bip] = 0; 
	}
    }
  
  for(auto &bipOccPair : bip2Occ)
    {
      auto &bip = std::get<0>(bipOccPair); 
      auto &occ = std::get<1>(bipOccPair); 
      for(auto bipHash : bipHashes)
	occ += bipHash.getPresence(bip).count(); 
    }

  // sort the bipartitions 
  auto sortBySecond = [](const std::pair<Bipartition,nat> &elemA, const std::pair<Bipartition,nat> &elemB)
    { 
      return elemA.second < elemB.second; 
    }; 

  auto bip2occVect = std::vector<std::pair<Bipartition,nat> >{bip2Occ.begin(), bip2Occ.end()}; 
  std::sort(bip2occVect.begin(), bip2occVect.end(), sortBySecond); 

  auto consensus = std::vector<Bipartition>{}; 
  auto minorityBips = std::vector<Bipartition>{}; 

  for(auto &bipPair : bip2occVect)
    {
      if(bipPair.first.count() == 1 )
      	continue; 
      else if(absThreshold <= bipPair.second )
	consensus.push_back(bipPair.first); 
      else 
	minorityBips.push_back(bipPair.first); 
    }

  if(isMRE)
    consensus = getRefinedConsensusTree(consensus, minorityBips); 

  assert(consensus.size() <= maxBip ); 

  result = bipEx.bipartitionsToTreeString(consensus, true); 

  return result; 
}


