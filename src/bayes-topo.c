#include "globals.h"
#include "common.h"

#include "config.h"
#include "axml.h"

#include "main-common.h"




topol  *setupTopol (int maxtips)
{
  topol   *tpl;

  if (! (tpl = (topol *) malloc(sizeof(topol))) || 
      ! (tpl->links = (connptr) malloc((2*maxtips-3) * sizeof(connect))))
    {
      printf("ERROR: Unable to get topology memory");
      tpl = (topol *) NULL;
    }
  else 
    {
      tpl->likelihood  = unlikely;
      tpl->start       = (node *) NULL;
      tpl->nextlink    = 0;
      tpl->ntips       = 0;
      tpl->nextnode    = 0;    
      tpl->scrNum      = 0;     /* position in sorted list of scores */
      tpl->tplNum      = 0;     /* position in sorted list of trees */	      
    }
  
  return  tpl;
}


void  freeTopol (topol *tpl)
{
  free(tpl->links);
  free(tpl);
}


static int  cmpTipVal (void *v1, void *v2)
{
  int  i1, i2;
  
  i1 = *((int *) v1);
  i2 = *((int *) v2);
  return  (i1 < i2) ? -1 : ((i1 == i2) ? 0 : 1);
}


static void  *tipValPtr (nodeptr p)
{ 
  return  (void *) & p->number;
}


static nodeptr minSubtreeTip (nodeptr  p0, int numsp)
{ 
  nodeptr  minTip, p, testTip;

  if (isTip(p0->number, numsp)) 
    return p0;

  p = p0->next;

  minTip = minSubtreeTip(p->back, numsp);

  while ((p = p->next) != p0) 
    {
      testTip = minSubtreeTip(p->back, numsp);
      if (cmpTipVal(tipValPtr(testTip), tipValPtr(minTip)) < 0)
        minTip = testTip;
    }
  return minTip;
} 

static nodeptr  minTreeTip (nodeptr  p, int numsp)
{
  nodeptr  minp, minpb;

  minp  = minSubtreeTip(p, numsp);
  minpb = minSubtreeTip(p->back, numsp);
  return (cmpTipVal(tipValPtr(minp), tipValPtr(minpb)) < 0 ? minp : minpb);
}

static int saveSubtree (nodeptr p, topol *tpl, int numsp, int numBranches)  
{
  connptr  r, r0;
  nodeptr  q, s;
  int      t, t0, t1, k;

  r0 = tpl->links;
  r = r0 + (tpl->nextlink)++;
  r->p = p;
  r->q = q = p->back;

  for(k = 0; k < numBranches; k++)
    r->z[k] = p->z[k];

  r->descend = 0;                     /* No children (yet) */

  if (isTip(q->number, numsp)) 
    {
      r->valptr = tipValPtr(q);         /* Assign value */
    }
  else 
    {                              /* Internal node, look at children */
      s = q->next;                      /* First child */
      do 
	{
	  t = saveSubtree(s, tpl, numsp, numBranches);        /* Generate child's subtree */

	  t0 = 0;                         /* Merge child into list */
	  t1 = r->descend;
	  while (t1 && (cmpTipVal(r0[t1].valptr, r0[t].valptr) < 0)) {
	    t0 = t1;
	    t1 = r0[t1].sibling;
          }
	  if (t0) r0[t0].sibling = t;  else  r->descend = t;
	  r0[t].sibling = t1;

	  s = s->next;                    /* Next child */
        } while (s != q);

      r->valptr = r0[r->descend].valptr;   /* Inherit first child's value */
      }                                 /* End of internal node processing */

  return  (r - r0);
}

void saveTree (tree *tr, topol *tpl)
/*  Save a tree topology in a standard order so that first branches
 *  from a node contain lower value tips than do second branches from
 *  the node.  The root tip should have the lowest value of all.
 */
{
  connptr  r;  
  
  tpl->nextlink = 0;                             /* Reset link pointer */
  r = tpl->links + saveSubtree(minTreeTip(tr->start, tr->mxtips), tpl, tr->mxtips, tr->numBranches);  /* Save tree */
  r->sibling = 0;
  
  tpl->likelihood = tr->likelihood;
  tpl->start      = tr->start;
  tpl->ntips      = tr->ntips;
  tpl->nextnode   = tr->nextnode;    
  
}



boolean restoreTree (topol *tpl, tree *tr)
{ 
  connptr  r;
  nodeptr  p, p0;    
  int  i;

  for (i = 1; i <= 2*(tr->mxtips) - 2; i++) 
    {  
      /* Uses p = p->next at tip */
      p0 = p = tr->nodep[i];
      do 
	{
	  p->back = (nodeptr) NULL;
	  p = p->next;
	} 
      while (p != p0);
    }

  /*  Copy connections from topology */

  for (r = tpl->links, i = 0; i < tpl->nextlink; r++, i++)     
    hookup(r->p, r->q, r->z, tr->numBranches);      

  tr->likelihood = tpl->likelihood;
  tr->start      = tpl->start;
  tr->ntips      = tpl->ntips;
  
  tr->nextnode   = tpl->nextnode;    

  evaluateGenericWrapper(tr, tr->start, TRUE);
  return TRUE;
}
