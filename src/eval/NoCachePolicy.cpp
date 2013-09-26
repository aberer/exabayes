#include "NoCachePolicy.hpp"
#include <cassert>

NoCachePolicy::NoCachePolicy(const TreeAln &traln )
{
}


std::unique_ptr<ArrayPolicy> NoCachePolicy::clone() const  
{
  return std::unique_ptr<ArrayPolicy>(new NoCachePolicy(*this) ); 
}


void NoCachePolicy::accountForRejectionPolicy(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat>& invalidNodes, ArrayOrientation &arrayOrient)
{
#if 1 
  // auto result = std::vector<nat>{}; 
  // for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
  //   if(partitions[i])
  //     result.push_back(i); 

  // tout << "accounting for rejection of " << invalidNodes
  //      << " for partitons " << result << std::endl; 
#endif

  for(nat i = 0; i < partitions.size( ); ++i)
    {
      if(not partitions[i])
	continue; 
      
      for(auto &elem : invalidNodes)
	{
	  nat id = elem - traln.getNumberOfTaxa()  -1 ; 
	  arrayOrient.setInvalid(i,id); 
	}
    }
}  
