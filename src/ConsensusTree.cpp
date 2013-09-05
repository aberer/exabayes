#include "ConsensusTree.hpp"
#include <map> 
#include <cassert> 
#include <algorithm>
 
ConsensusTree::ConsensusTree(std::string file )
  : bipEx(std::vector<std::string>{file})
{
  bipEx.extractBipsNew(); 
  const auto &hash = bipEx.getBipartitionHashes()[0]; 
  totalTrees = hash.getTreesAdded();
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
  assert(bipHashes.size() == 1 ); 
  const auto& bipHash = bipHashes[0]; 
  auto bip2Occ = std::vector<std::pair<Bipartition,nat> >{}; 
  nat maxBip = bipEx.getTaxa().size() - 3; 

  nat absThreshold = 0; 
  if(isMRE)
    absThreshold = nat(bipHash.getTreesAdded() * .5);
  else 
    absThreshold = nat(bipHash.getTreesAdded() * threshold); 

  assert(absThreshold > 0. ); 

  for(auto bip : bipHash)
    {
      auto cpy = bip.first ; 
      auto pair = std::make_pair(cpy, bipHash.getPresence(cpy).count()); 
      bip2Occ.push_back(pair);
    }
  
  // sort the bipartitions 
  auto sortBySecond = [](const std::pair<Bipartition,nat> &elemA, const std::pair<Bipartition,nat> &elemB)
    { 
      return elemA.second < elemB.second; 
    }; 
  std::sort(bip2Occ.begin(), bip2Occ.end(), sortBySecond); 

  auto consensus = std::vector<Bipartition>{}; 
  auto minorityBips = std::vector<Bipartition>{}; 

  for(auto &bipPair : bip2Occ)
    {
      if(bipPair.first.count() == 1 )
	continue; 
      else if(absThreshold < bipPair.second )
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


