#include <cassert>

#include "Topology.hpp"

// deprecated 

Topology::Topology(int numTax)  
{  
  std::cout << "initializing deprecated topology  " << std::endl; 
}


void Topology::saveTopology(TreeAln &traln)
{
  branches.clear(); 
  branches = traln.extractBranches();
}


void Topology::restoreTopology(TreeAln &traln)
{
  traln.unlinkTree();
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 

  for(auto b : branches)
    {
      nodeptr p = traln.getUnhookedNode(b.getPrimNode()) ,
	q = traln.getUnhookedNode(b.getSecNode()); 
      
      double tmp = b.getLength(); 
      traln.clipNode(p,q,tmp ); 
      
    }  
}

