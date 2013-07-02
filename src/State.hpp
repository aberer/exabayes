/**
   @brief represents a state in the space we are integrating over. 

   This compises all relevant parameters and the topology. 
 */ 


#ifndef __STATE_H
#define __STATE_H

#include "GlobalVariables.hpp"
#include "Partition.hpp"
#include "Topology.hpp"


class State
{
public: 
  State(const TreeAln &traln)
    : topology(traln.getTr()->mxtips)    
  {
    nat numPart = traln.getNumberOfPartitions(); 
    for(nat i = 0; i < numPart; ++i)
      partitions.push_back(Partition(i,traln));     
  }
  
  Topology& accessTopology(){ return topology; }
  Partition& accessPartition(int num) { return partitions[num]; }

private: 
  Topology topology; 
  vector<Partition> partitions; 
}; 


#endif
