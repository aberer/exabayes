#include "TreeRandomizer.hpp"
#include "TreeAln.hpp"
#include "treeRead.h"


TreeRandomizer::TreeRandomizer(int seed)
  : rand(seed)
{
}


void TreeRandomizer::randomizeTree(TreeAlnPtr &traln)
{
  tree *tr = traln->getTr();

  nodeptr 
    p, 
    f, 
    randomBranch,
    *branches = (nodeptr *)exa_malloc(sizeof(nodeptr) * (2 * tr->mxtips));    
  
  int 
    nextsp, 
    *perm = (int *) exa_malloc((tr->mxtips + 1) * sizeof(int)), 
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
  
  buildSimpleTreeRandom(traln, perm[1], perm[2], perm[3]);
  
  while(tr->ntips < tr->mxtips) 
    {	             
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];            
      
      buildNewTip(traln,p);  	
      
      f = findAnyTip(tr->start, tr->mxtips);
      f = f->back;
      
      branchCounter = 1;
      branches[0] = f;
      markBranches(traln, branches, f, &branchCounter, tr->mxtips);

      assert(branchCounter == ((2 * (tr->ntips - 1)) - 3));
      
      randomBranch = branches[rand.drawRandInt(branchCounter)];
      
      insertTaxon( traln, p->back, randomBranch);      
    }
  
  exa_free(perm); 
  exa_free(branches);  
}


int TreeRandomizer::markBranches(TreeAlnPtr &traln, nodeptr *branches, nodeptr p, int *counter, int numsp)
{
  if(isTip(p->number, numsp))
    return 0;
  else
    {
      branches[*counter] = p->next;
      branches[*counter + 1] = p->next->next;
      
      *counter = *counter + 2;
      
      return ((2 + markBranches(traln, branches, p->next->back, counter, numsp) + 
	       markBranches(traln, branches, p->next->next->back, counter, numsp)));
    }
}


void TreeRandomizer::insertTaxon (TreeAlnPtr &traln, nodeptr p, nodeptr q)
{
  nodeptr  r;
  r = q->back;
  
  traln->clipNodeDefault(  p->next,       q);
  traln->clipNodeDefault( p->next->next, r); 
}


nodeptr TreeRandomizer::buildNewTip (TreeAlnPtr &traln, nodeptr p)
{ 
  nodeptr  q;
  tree *tr = traln->getTr();

  q = tr->nodep[(tr->nextnode)++];
  traln->clipNodeDefault( p, q);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 


void TreeRandomizer::buildSimpleTreeRandom (TreeAlnPtr &traln, int ip, int iq, int ir)
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
  
  traln->clipNodeDefault( p, tr->nodep[iq]);
  
  s = buildNewTip( traln, tr->nodep[ir]);
  
  insertTaxon(traln, s, p);
}

