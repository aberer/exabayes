

#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "adapters.h"



/* RANDOM TREE STUFF  */

static void insertTaxon (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;

  r = q->back;
  
  exa_hookupDefault(tr, p->next,       q);
  exa_hookupDefault(tr, p->next->next, r); 
} 

static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  exa_hookupDefault(tr, p, q);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTreeRandom (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  
    p, 
    s;
  
  int  
    i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  
  exa_hookupDefault(tr, p, tr->nodep[iq]);
  
  s = buildNewTip(tr, tr->nodep[ir]);
  
  insertTaxon(tr, s, p);
}


static void makePermutation(int *perm, int n, tree *tr)
{    
  int  
    i, 
    j, 
    k;    

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {    
      k =  drawGlobalRandIntBound(n + 1 - i);

      assert(i + k <= n);
      
      j        = perm[i];
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

static int markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp)
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



void exa_makeRandomTree(tree *tr)
{  
  nodeptr 
    p, 
    f, 
    randomBranch,
    *branches = (nodeptr *)malloc(sizeof(nodeptr) * (2 * tr->mxtips));    
  
  int 
    nextsp, 
    *perm = (int *)malloc((tr->mxtips + 1) * sizeof(int)), 
    branchCounter;                      
  
  makePermutation(perm, tr->mxtips, tr);              
  
  tr->ntips = 0;       	       
  tr->nextnode = tr->mxtips + 1;    
  
  buildSimpleTreeRandom(tr, perm[1], perm[2], perm[3]);
  
  while(tr->ntips < tr->mxtips) 
    {	             
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];            
      
      buildNewTip(tr, p);  	
      
      f = findAnyTip(tr->start, tr->mxtips);
      f = f->back;
      
      branchCounter = 1;
      branches[0] = f;
      markBranches(branches, f, &branchCounter, tr->mxtips);

      assert(branchCounter == ((2 * (tr->ntips - 1)) - 3));
      
      randomBranch = branches[drawGlobalRandIntBound(branchCounter)];
      
      insertTaxon(tr, p->back, randomBranch);      
    }
  
  exa_free(perm); 
  exa_free(branches);
}

/* END RANDOM TRE  */

