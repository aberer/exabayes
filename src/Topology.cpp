
#include "Topology.hpp"
#include "output.h"


// #define TOPO_INFO  


Topology::Topology(int numTax)  
{  
  for(int i = 0  ; i < 2  * numTax - 3; ++i)
    {
      branch *b =  (branch*)exa_calloc(1,sizeof(branch)); 
      branches.push_back(b);
    }
}


Topology::~Topology()
{
  for(auto b : branches)
    exa_free(b);  
}


void Topology::traverseAndSave(TreeAln &traln, nodeptr p, nat &number)
{
  assert(traln.getNumBranches() == 1); 
  tree *tr = traln.getTr();

#ifdef TOPO_INFO
  cout << "saving branch " << p->number << " and " << p->back->number << " (" << traln.getBranchLength(p->number,0) << ") " << endl; 
#endif

  // save branch 
  branch *b = branches[number]; 
  b->thisNode = p->number; 
  b->thatNode = p->back->number; ; 
  b->length[0] = traln.getBranchLength( p,0);   
  number++; 

  if(NOT isTip(p->number, tr->mxtips))
    {
      traverseAndSave(traln, p->next->back, number); 
      traverseAndSave(traln, p->next->next->back, number);   
    }
}




void Topology::saveTopology(TreeAln &traln)
{
  tree *tr = traln.getTr(); 

  // clear all saved info 
  for(branch* b : this->branches )
    {
      b->thisNode = 0; 
      b->thatNode = 0; 
      b->length[0] = 0; 
    }

  assert(traln.getNumBranches()== 1 );
  nat ctr = 0; 
  traverseAndSave(traln, tr->start->back,ctr);  
  assert(ctr == branches.size()); 
}


void Topology::restoreTopology(TreeAln &traln)
{
  traln.unlinkTree();
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 

  for(branch *b : this->branches)
    {
#ifdef TOPO_INFO
      cout << "hooking " << b->thisNode << " / " << b->thatNode << endl; 
#endif
      nodeptr p = traln.getUnhookedNode(b->thisNode),
	q = traln.getUnhookedNode(b->thatNode); 
      
      hookup(p, q, b->length, numBranches); 
    }

  debug_checkTreeConsistency(traln.getTr());
}

