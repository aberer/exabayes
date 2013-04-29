#include "TreeRandomizer.hpp"
#include "TreeAln.hpp"
#include "adapters.h"


TreeRandomizer::TreeRandomizer(int seed, TreeAln *_traln)
  : rand(seed)
  , traln(_traln)
{
}


void TreeRandomizer::randomizeTree()
{
  tree *tr = traln->getTr();

  nodeptr 
    p, 
    f, 
    randomBranch,
    *branches = (nodeptr *)malloc(sizeof(nodeptr) * (2 * tr->mxtips));    
  
  int 
    nextsp, 
    *perm = (int *)malloc((tr->mxtips + 1) * sizeof(int)), 
    branchCounter;                      

  // permutation 
  {
  for (int i = 1; i <= tr->mxtips; i++)
    perm[i] = i;

  
  for (int i = 1; i <= tr->mxtips; i++)
    {
      int k = rand.drawRandInt(tr->mxtips + 1 -i);

      assert(i + k <= tr->mxtips);

      int j = perm[i];
      perm[i] = perm[i + k];
      perm[i + k] = j;
    }
  }

  
  tr->ntips = 0;       	       
  tr->nextnode = tr->mxtips + 1;    
  
  buildSimpleTreeRandom(perm[1], perm[2], perm[3]);
  
  while(tr->ntips < tr->mxtips) 
    {	             
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];            
      
      buildNewTip(p);  	
      
      f = findAnyTip(tr->start, tr->mxtips);
      f = f->back;
      
      branchCounter = 1;
      branches[0] = f;
      markBranches(branches, f, &branchCounter, tr->mxtips);

      assert(branchCounter == ((2 * (tr->ntips - 1)) - 3));
      
      randomBranch = branches[rand.drawRandInt(branchCounter)];
      
      insertTaxon( p->back, randomBranch);      
    }
  
  exa_free(perm); 
  exa_free(branches);
  
}


int TreeRandomizer::markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp)
{
  if(isTip(p->number, numsp))
    return 0;
  else
    {
      branches[*counter] = p->next;
      branches[*counter + 1] = p->next->next;
      
      *counter = *counter + 2;
      
      return ((2 + markBranches(branches, p->next->back, counter, numsp) + 
	       markBranches(branches, p->next->next->back, counter, numsp)));
    }
}


void TreeRandomizer::insertTaxon (nodeptr p, nodeptr q)
{
  nodeptr  r;

  r = q->back;
  tree *tr = traln->getTr();
  
  exa_hookupDefault(tr,  p->next,       q);
  exa_hookupDefault(tr, p->next->next, r); 
}


nodeptr TreeRandomizer::buildNewTip (nodeptr p)
{ 
  nodeptr  q;
  tree *tr = traln->getTr();

  q = tr->nodep[(tr->nextnode)++];
  exa_hookupDefault(tr, p, q);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 


void TreeRandomizer::buildSimpleTreeRandom (int ip, int iq, int ir)
{    
  nodeptr  
    p, 
    s;
  
  int  
    i;

  tree *tr = traln->getTr();
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  
  exa_hookupDefault(tr, p, tr->nodep[iq]);
  
  s = buildNewTip( tr->nodep[ir]);
  
  insertTaxon(s, p);
}

