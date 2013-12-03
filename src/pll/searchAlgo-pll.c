/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include "mem_alloc.h"

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>



#include "axml.h"


/** @file searchAlgo.c
    
    @brief Collection of routines for performing likelihood computation and branch optimization.

    Detailed description to appear soon.
*/




extern double accumulatedTime;

extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern char run_id[128];
extern double masterTime;
extern partitionLengths pLengths[MAX_MODEL];
extern char binaryCheckpointName[1024];
extern char binaryCheckpointInputName[1024];

boolean initrav (tree *tr, partitionList *pr, nodeptr p)
{ 
  nodeptr  q;

  if (!isTip(p->number, tr->mxtips)) 
  {      
    q = p->next;

    do 
    {	   
      if (! initrav(tr, pr, q->back))  return FALSE;
      q = q->next;	
    } 
    while (q != p);  

    newviewGeneric(tr, pr, p, FALSE, NULL);
  }

  return TRUE;
} 





/** @brief Optimize the length of a specific branch

    Optimize the length of the branch connecting \a p and \a p->back
    for each partition (\a tr->numBranches) in tree \a tr.
 
    @param tr
      The tree structure

    @param pr
      Partition list
 
    @param p
      Endpoints of branch to be optimized 
*/
void update(tree *tr, partitionList *pr, nodeptr p)
{       
  nodeptr  q; 
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
  q = p->back;   

  for(i = 0; i < numBranches; i++)
    z0[i] = q->z[i];    

  double fake;  
  double fakeL = 0; 
  double fakeM = 0;  

  if(numBranches > 1)
    makenewzGeneric(tr, pr, p, q, z0, newzpercycle, z, &fakeM, &fake, fakeL , TRUE, NULL);
  else
    makenewzGeneric(tr, pr, p, q, z0, newzpercycle, z, &fakeM, &fake, fakeL, FALSE, NULL);

  for(i = 0; i < numBranches; i++)
  {         
    if(!tr->partitionConverged[i])
    {	  
      if(ABS(z[i] - z0[i]) > deltaz)  
      {	      
        tr->partitionSmoothed[i] = FALSE;
      }	 

      p->z[i] = q->z[i] = z[i];	 
    }
  }
}

/** @brief Branch length optimization of specific branches

    Optimize the length of branches that have \a p as an endpoint 

    @param tr
      The tree structure

    @param pr
      Partition list

    @param p
      Endpoint of branches to be optimized
*/
void smooth (tree *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  update(tr, pr, p);    /*  Adjust branch */

  if (! isTip(p->number, tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      smooth(tr, pr, q->back);
      q = q->next;
    }	

    if(numBranches > 1 && !tr->useRecom)
      newviewGeneric(tr, pr,p, TRUE,NULL);
    else
      newviewGeneric(tr, pr,p, FALSE, NULL);
  }
} 

/**  @brief Check whether the branches in all partitions have been optimized
 
     Check if all branches in all partitions have reached the threshold for
     optimization. If at least one branch can be optimized further return \b FALSE.

     @param tr
       The tree structure

     @return
       If at least one branch can be further optimized return \b FALSE,
       otherwise \b TRUE.
             
*/
static boolean allSmoothed(tree *tr, int numBranches)
{
  int i;
  boolean result = TRUE;

  for(i = 0; i < numBranches; i++)
  {
    if(tr->partitionSmoothed[i] == FALSE)
      result = FALSE;
    else
      tr->partitionConverged[i] = TRUE;
  }

  return result;
}


/** @brief Wrapper function for branch length optimization of the tree
  
    Perform \a maxtimes rounds of branch length optimization by running smooth()
    on all neighbour nodes of node \a tr->start.

    @param tr
      The tree structure

    @param maxtimes
      Number of optimization rounds to perform
*/
/* do maxtimes rounds of branch length optimization */
void smoothTree (tree *tr, partitionList *pr, int maxtimes)
{
	nodeptr  p, q;
	int i, count = 0;
    int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

	p = tr->start;
	for(i = 0; i < numBranches; i++)
		tr->partitionConverged[i] = FALSE;

	while (--maxtimes >= 0)
	{
		for(i = 0; i < numBranches; i++)
			tr->partitionSmoothed[i] = TRUE;

		smooth(tr, pr, p->back);
		if (!isTip(p->number, tr->mxtips))
		{
			q = p->next;
			while (q != p)
			{
				smooth(tr, pr, q->back);
				q = q->next;
			}
		}
		count++;

		if (allSmoothed(tr, numBranches)) break;
	}

	for(i = 0; i < numBranches; i++)
		tr->partitionConverged[i] = FALSE;
} 

/** @brief Optimize the branch length of edges around a specific node
    
    Optimize \a maxtimes the branch length of all (3) edges around a given node 
    \a p of a tree \a tr.

    @param tr
      The tree structure

    @param p
      The node around which to optimize the edges

    @param maxtimes
      Number of optimization rounds to perform
*/
void localSmooth (tree *tr, partitionList *pr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  int i;
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
  if (isTip(p->number, tr->mxtips)) return;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionConverged[i] = FALSE;	

  while (--maxtimes >= 0) 
  {     
    for(i = 0; i < NUM_BRANCHES; i++)
      tr->partitionSmoothed[i] = TRUE;

    q = p;
    do 
    {
      update(tr, pr, q);
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr, numBranches))
      break;
  }

  for(i = 0; i < NUM_BRANCHES; i++)
  {
    tr->partitionSmoothed[i] = FALSE; 
    tr->partitionConverged[i] = FALSE;
  }
}





/** @brief Reset an \a infoList

    Resets an \a infoList by setting elements \a node and \a likelihood
    of each element of the \a bestInfo list structure to \b NULL and
    \a unlikely, respectively.

    @param iList
      The given \a infoList.
*/
static void resetInfoList(infoList *iList)
{
  int 
    i;

  iList->valid = 0;

  for(i = 0; i < iList->n; i++)    
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = unlikely;
  }    
}

/** @brief Initialize an \a infoList

    Initialize an \a infoList by creating a \a bestInfo list structure
    of \a n elements and setting the attributes \a node and \a likelihood
    of each element of the \a bestInfo list structure to \b NULL and
    \a unlikely, respectively.

    @param iList
      The given \a infoList.

    @param n
      Number of elements to be created in the \a bestInfo list.
*/
static void initInfoList(infoList *iList, int n)
{
  int 
    i;

  iList->n = n;
  iList->valid = 0;
  iList->list = (bestInfo *)rax_malloc(sizeof(bestInfo) * (size_t)n);

  for(i = 0; i < n; i++)
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = unlikely;
  }
}

/** @brief Deallocate the contents of an \a infoList
    
    Deallocate the contents of a given \a infoList by freeing
    the memory used by its \a bestInfo list structure.

    @param iList
      The \a infoList to be used.
*/
static void freeInfoList(infoList *iList)
{ 
  rax_free(iList->list);   
}


/** @brief Insert a record in an \a infoList

    Insert the pair \a likelihood and \node into the \a bestList \a list
    of \a iList \b only if there already exists a pair in \a list 
    which has the \a likelihood attribute smaller than the given \a 
    likelihoodby. The insertion is done by replacing the smallest
    likelihood pair with the new pair.

    @param node
      The given node

    @param likelihood
      The given likelihood

    @param iList
      The given \a infoList where the record will possibly be appended.
*/
static void insertInfoList(nodeptr node, double likelihood, infoList *iList)
{
  int 
    i,
    min = 0;

  double 
    min_l =  iList->list[0].likelihood;

  for(i = 1; i < iList->n; i++)
  {
    if(iList->list[i].likelihood < min_l)
    {
      min = i;
      min_l = iList->list[i].likelihood;
    }
  }

  if(likelihood > min_l)
  {
    iList->list[min].likelihood = likelihood;
    iList->list[min].node = node;
    iList->valid += 1;
  }

  if(iList->valid > iList->n)
    iList->valid = iList->n;
}


/** @brief  Optimize branch lengths of region

    Optimize the branch lenghts of only a specific region. The branch optimization starts
    at a node \a p and is carried out in all nodes with distance upto \a region from \a p.

    @param tr
      The tree structure.
    
    @param p
      Node to start branch optimization from.

    @param region
      The allowed node distance from \p for which to still perform branch optimization.
*/
void smoothRegion (tree *tr, partitionList *pr, nodeptr p, int region)
{ 
  nodeptr  q;

  update(tr, pr, p);   /* Adjust branch */

  if (region > 0)
  {
    if (!isTip(p->number, tr->mxtips)) 
    {                                 
      q = p->next;
      while (q != p) 
      {
        smoothRegion(tr, pr, q->back, --region);
        q = q->next;
      }	

      newviewGeneric(tr, pr,p, FALSE, NULL);
    }
  }
}

/** @brief Wrapper function for optimizing the branch length of a region \a maxtimes times

    Optimize the branch lengths of a specific region \a maxtimes times. The branch optimization
    starts at a given node \a p and is carried out in all nodes with distance upto \a region
    from \a p.

    @param tr
      The tree structure.

    @param p
      Node to start branch optimization from.

    @param maxtimes
      Number of times to perform branch optimization.

    @pram region
      The allwed node distance from \p for which to still perform branch optimization.

    @todo
      In the previous version (before the model-sep merge) the loops were controlled by tr->numBranches,
      and now they are controlled by a constant NUM_BRANCHES. What is right?
*/
void regionalSmooth (tree *tr, partitionList *pr, nodeptr p, int maxtimes, int region)
{
  nodeptr  q;
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if (isTip(p->number, tr->mxtips)) return;            /* Should be an error */

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionConverged[i] = FALSE;

  while (--maxtimes >= 0) 
  {	
    for(i = 0; i < NUM_BRANCHES; i++)
      tr->partitionSmoothed[i] = TRUE;

    q = p;
    do 
    {
      smoothRegion(tr, pr, q, region);
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr, numBranches))
      break;
  }

  for(i = 0; i < NUM_BRANCHES; i++) {
    tr->partitionSmoothed[i] = FALSE;
    tr->partitionConverged[i] = FALSE;
  }
} 

/* @brief Split the tree into two components and optimize new branch length

   Split the tree into two components. The disconnection point is node \a p.
   First, a branch length is computed for the newly created branch between nodes
   \a p->next->back and \a p->next->next->back and then the two nodes are
   connected (hookup). Disconnection is done by setting \a p->next->next->back
   and \a p->next->back to \b NULL.

   @param tr
     The tree structure

   @param p
     The node at which the tree should be decomposed into two components.

   @param numBranches
     Number of branches per partition

   @return q
     Node from the disconnected component

   @todo
     Why do we return this node?
*/
nodeptr  removeNodeBIG (tree *tr, partitionList *pr, nodeptr p, int numBranches)
{  
  double   zqr[numBranches], result[numBranches];
  nodeptr  q, r;
  int i;

  q = p->next->back;
  r = p->next->next->back;

  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        

  double fake =0; 
  double fakeM = 0; 
  double fakeL = 0; 
  makenewzGeneric(tr, pr, q, r, zqr, iterations, result, &fakeM, &fake, fakeL , FALSE, NULL);

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 

  p->next->next->back = p->next->back = (node *) NULL;

  return  q; 
}


/** @brief Split the tree into two components and recompute likelihood

    Split the tree into two component. The disconnection point is node \a p.
    Set the branch length of the new node between \a p->next->back and
    \a p->next->next->back to \a tr->currentZQR and then decompose the tree
    into two components by setting \a p->next->back and \a p->next->next->back
    to \b NULL.

    @param tr
      The tree structure

    @param p
      The node at which the tree should be decomposed into two components.

    @return q
      the node after \a p

    @todo
      Why do we return this node? Why do we set to tr->currentZQR and not compute
      new optimized length? What is tr->currentZQR? 
*/
nodeptr  removeNodeRestoreBIG (tree *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q, r;

  q = p->next->back;
  r = p->next->next->back;  

  newviewGeneric(tr, pr,q, FALSE, NULL);
  newviewGeneric(tr, pr,r, FALSE, NULL);

  hookup(q, r, tr->currentZQR, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  p->next->next->back = p->next->back = (node *) NULL;

  return  q;
}



/* @brief

   @todo
     What is tr->lzi ? What is thorough insertion?
*/
boolean insertBIG (tree *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  r = q->back;
  s = p->back;

  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];

  if(tr->thoroughInsertion)
  { 
    double  zqr[numBranches], zqs[numBranches], zrs[numBranches], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;
    double defaultArray[numBranches];
    double e1[numBranches], e2[numBranches], e3[numBranches];
    double *qz;

    qz = q->z;

    for(i = 0; i < numBranches; i++)
      defaultArray[i] = defaultz;

    double fakeL = 0; 
    double fake = 0; 
    double fakeM = 0; 

    makenewzGeneric(tr, pr, q, r, qz, iterations, zqr,  &fakeM, &fake , fakeL, FALSE, NULL);
    /* the branch lengths values will be estimated using q, r and s
     * q-s are not connected, but both q and s have a valid LH vector , so we can call makenewzGeneric  to get a value for
     * lzsum, which is then use to generate reasonable starting values e1, e2, e3 for the new branches we create after the       insertion
     */


    /* double fakeM = 0;  */

    makenewzGeneric(tr, pr, q, s, defaultArray, iterations, zqs, &fakeM , &fake, fakeL , FALSE, NULL);
    makenewzGeneric(tr, pr, r, s, defaultArray, iterations, zrs, &fakeM, &fake, fakeL , FALSE, NULL);


    for(i = 0; i < numBranches; i++)
    {
      lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
      lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
      lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
      lzsum = 0.5 * (lzqr + lzqs + lzrs);

      lzq = lzsum - lzrs;
      lzr = lzsum - lzqs;
      lzs = lzsum - lzqr;
      lzmax = log(zmax);

      if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
      else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
      else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          

      e1[i] = exp(lzq);
      e2[i] = exp(lzr);
      e3[i] = exp(lzs);
    }
    hookup(p->next,       q, e1, numBranches);
    hookup(p->next->next, r, e2, numBranches);
    hookup(p,             s, e3, numBranches);      		  
  }
  else
  {       
    double  z[numBranches];

    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      

      if(z[i] < zmin) 
        z[i] = zmin;
      if(z[i] > zmax)
        z[i] = zmax;
    }

    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }

  newviewGeneric(tr, pr,p, FALSE, NULL);

  if(tr->thoroughInsertion)
  {     
    localSmooth(tr, pr, p, MAX_LOCAL_SMOOTHING_ITERATIONS);
    for(i = 0; i < numBranches; i++)
    {
      tr->lzq[i] = p->next->z[i];
      tr->lzr[i] = p->next->next->z[i];
      tr->lzs[i] = p->z[i];            
    }
  }           
  return  TRUE;
}


boolean insertRestoreBIG (tree *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;

  r = q->back;
  s = p->back;

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, numBranches);
    hookup(p->next->next, r, tr->currentLZR, numBranches);
    hookup(p,             s, tr->currentLZS, numBranches);
  }
  else
  {       
    double  z[NUM_BRANCHES];
    int i;

    for(i = 0; i < numBranches; i++)
    {
      double zz;
      zz = sqrt(q->z[i]);     
      if(zz < zmin) 
        zz = zmin;
      if(zz > zmax)
        zz = zmax;
      z[i] = zz;
    }

    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }   

  newviewGeneric(tr, pr,p, FALSE, NULL);

  return  TRUE;
}


static void restoreTopologyOnly(tree *tr, bestlist *bt, int numBranches)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[NUM_BRANCHES], pz[NUM_BRANCHES], p1z[NUM_BRANCHES], p2z[NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;

  p1 = p->next->back;
  p2 = p->next->next->back;

  //memcpy(p1z, p1->z, numBranches*sizeof(double));
  //memcpy(p2z, p2->z, numBranches*sizeof(double));
  //memcpy(qz, q->z, numBranches*sizeof(double));
  //memcpy(pz, p->z, numBranches*sizeof(double));
  for(i = 0; i < numBranches; i++)
  {
    p1z[i] = p1->z[i];
    p2z[i] = p2->z[i];
  }

  hookup(p1, p2, tr->currentZQR, numBranches);

  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }

  r = q->back;
  s = p->back;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, numBranches);
    hookup(p->next->next, r, tr->currentLZR, numBranches);
    hookup(p,             s, tr->currentLZS, numBranches);
  }
  else
  { 	
    double  z[NUM_BRANCHES];	
    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      
      if(z[i] < zmin)
        z[i] = zmin;
      if(z[i] > zmax)
        z[i] = zmax;
    }
    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }     

  tr->likelihood = tr->bestOfNode;

  saveBestTree(bt, tr, numBranches);

  tr->likelihood = currentLH;

  hookup(q, r, qz, numBranches);

  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(tr->thoroughInsertion)    
    hookup(p, s, pz, numBranches);

  hookup(p->next,       p1, p1z, numBranches);
  hookup(p->next->next, p2, p2z, numBranches);
}


boolean testInsertBIG (tree *tr, partitionList *pr, nodeptr p, nodeptr q)
{

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = TRUE;
  double startLH = tr->endLH;
  int i;

  r = q->back; 
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }



  if(doIt)
  {     
    if (! insertBIG(tr, pr, p, q))       return FALSE;

    evaluateGeneric(tr, pr, p->next->next, FALSE, NULL);

    if(tr->likelihood > tr->bestOfNode)
    {
      tr->bestOfNode = tr->likelihood;
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < numBranches; i++)
      {
        tr->currentZQR[i] = tr->zqr[i];           
        tr->currentLZR[i] = tr->lzr[i];
        tr->currentLZQ[i] = tr->lzq[i];
        tr->currentLZS[i] = tr->lzs[i];      
      }
    }

    if(tr->likelihood > tr->endLH)
    {			  
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < numBranches; i++)
        tr->currentZQR[i] = tr->zqr[i];      
      tr->endLH = tr->likelihood;                      
    }        

    hookup(q, r, qz, numBranches);

    p->next->next->back = p->next->back = (nodeptr) NULL;

    if(tr->thoroughInsertion)
    {
      nodeptr s = p->back;
      hookup(p, s, pz, numBranches);
    } 

    if((tr->doCutoff) && (tr->likelihood < startLH))
    {
      tr->lhAVG += (startLH - tr->likelihood);
      tr->lhDEC++;
      if((startLH - tr->likelihood) >= tr->lhCutoff)
        return FALSE;	    
      else
        return TRUE;
    }
    else
      return TRUE;
  }
  else
    return TRUE;  
}




void addTraverseBIG(tree *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
  {              
    if (! testInsertBIG(tr, pr, p, q))  return;

  }

  if ((!isTip(q->number, tr->mxtips)) && (--maxtrav > 0)) 
  {    
    addTraverseBIG(tr, pr, p, q->next->back, mintrav, maxtrav);
    addTraverseBIG(tr, pr, p, q->next->next->back, mintrav, maxtrav);
  }
} 





int rearrangeBIG(tree *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav)
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  boolean doP = TRUE, doQ = TRUE;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  q = p->back;




  if (!isTip(p->number, tr->mxtips) && doP) 
  {     
    p1 = p->next->back;
    p2 = p->next->next->back;


    if(!isTip(p1->number, tr->mxtips) || !isTip(p2->number, tr->mxtips))
    {
      for(i = 0; i < numBranches; i++)
      {
        p1z[i] = p1->z[i];
        p2z[i] = p2->z[i];	   	   
      }

      if (! removeNodeBIG(tr, pr, p,  numBranches)) return badRear;

      if (!isTip(p1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, p, p1->next->back,
            mintrav, maxtrav);         

        addTraverseBIG(tr, pr, p, p1->next->next->back,
            mintrav, maxtrav);          
      }

      if (!isTip(p2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, p, p2->next->back,
            mintrav, maxtrav);
        addTraverseBIG(tr, pr, p, p2->next->next->back,
            mintrav, maxtrav);          
      }

      hookup(p->next,       p1, p1z, numBranches);
      hookup(p->next->next, p2, p2z, numBranches);
      newviewGeneric(tr, pr,p, FALSE, NULL);
    }
  }  

  if (!isTip(q->number, tr->mxtips) && maxtrav > 0 && doQ) 
  {
    q1 = q->next->back;
    q2 = q->next->next->back;

    /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
      ((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
    if (
        (
         ! isTip(q1->number, tr->mxtips) && 
         (! isTip(q1->next->back->number, tr->mxtips) || ! isTip(q1->next->next->back->number, tr->mxtips))
        )
        ||
        (
         ! isTip(q2->number, tr->mxtips) && 
         (! isTip(q2->next->back->number, tr->mxtips) || ! isTip(q2->next->next->back->number, tr->mxtips))
        )
       )
    {

      for(i = 0; i < numBranches; i++)
      {
        q1z[i] = q1->z[i];
        q2z[i] = q2->z[i];
      }

      if (! removeNodeBIG(tr, pr, q, numBranches)) return badRear;

      mintrav2 = mintrav > 2 ? mintrav : 2;

      if (/*! q1->tip*/ !isTip(q1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, q, q1->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, pr, q, q1->next->next->back,
            mintrav2 , maxtrav);         
      }

      if (/*! q2->tip*/ ! isTip(q2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, q, q2->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, pr, q, q2->next->next->back,
            mintrav2 , maxtrav);          
      }	   

      hookup(q->next,       q1, q1z, numBranches);
      hookup(q->next->next, q2, q2z, numBranches);

      newviewGeneric(tr, pr,q, FALSE, NULL);
    }
  } 

  return  1;
} 





static double treeOptimizeRapid(tree *tr, partitionList *pr, int mintrav, int maxtrav, analdef *adef, bestlist *bt, infoList *iList)
{
  int i, index,
      *perm = (int*)NULL;   

  nodeRectifier(tr);



  if (maxtrav > tr->mxtips - 3)  
    maxtrav = tr->mxtips - 3;  



  resetInfoList(iList);

  resetBestTree(bt);

  tr->startLH = tr->endLH = tr->likelihood;

  if(tr->doCutoff)
  {
    if(tr->bigCutoff)
    {	  
      if(tr->itCount == 0)    
        tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
      else    		 
        tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
    }
    else
    {
      if(tr->itCount == 0)    
        tr->lhCutoff = tr->likelihood / -1000.0;    
      else    		 
        tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
    }    

    tr->itCount = tr->itCount + 1;
    tr->lhAVG = 0;
    tr->lhDEC = 0;
  }

  /*
     printf("DoCutoff: %d\n", tr->doCutoff);
     printf("%d %f %f %f\n", tr->itCount, tr->lhAVG, tr->lhDEC, tr->lhCutoff);

     printf("%d %d\n", mintrav, maxtrav);
     */

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
  {           
    tr->bestOfNode = unlikely;          

    if(adef->permuteTreeoptimize)
      index = perm[i];
    else
      index = i;     

    if(rearrangeBIG(tr, pr, tr->nodep[index], mintrav, maxtrav))
    {    
      if(tr->thoroughInsertion)
      {
        if(tr->endLH > tr->startLH)                 	
        {			   	     
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }
        else
        { 		  
          if(tr->bestOfNode != unlikely)
            restoreTopologyOnly(tr, bt, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }	   
      }
      else
      {
        insertInfoList(tr->nodep[index], tr->bestOfNode, iList);	    
        if(tr->endLH > tr->startLH)                 	
        {		      
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
        }	    	  
      }
    }     
  }     

  if(!tr->thoroughInsertion)
  {           
    tr->thoroughInsertion = TRUE;  

    for(i = 0; i < iList->valid; i++)
    { 	  
      tr->bestOfNode = unlikely;

      if(rearrangeBIG(tr, pr, iList->list[i].node, mintrav, maxtrav))
      {	  
        if(tr->endLH > tr->startLH)                 	
        {	 	     
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }
        else
        { 

          if(tr->bestOfNode != unlikely)
          {	     
            restoreTopologyOnly(tr, bt, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
          }	
        }      
      }
    }       

    tr->thoroughInsertion = FALSE;
  }

  if(adef->permuteTreeoptimize)
    rax_free(perm);

  return tr->startLH;     
}




boolean testInsertRestoreBIG (tree *tr, partitionList *pr, nodeptr p, nodeptr q)
{    
  if(tr->thoroughInsertion)
  {
    if (! insertBIG(tr, pr, p, q))       return FALSE;

    evaluateGeneric(tr, pr, p->next->next, FALSE, NULL);
  }
  else
  {
    if (! insertRestoreBIG(tr, pr, p, q))       return FALSE;

    {
      nodeptr x, y;
      x = p->next->next;
      y = p->back;

      if(! isTip(x->number, tr->mxtips) && isTip(y->number, tr->mxtips))
      {
        while ((! x->x)) 
        {
          if (! (x->x))
            newviewGeneric(tr, pr,x, FALSE, NULL);
        }
      }

      if(isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! y->x)) 
        {		  
          if (! (y->x))
            newviewGeneric(tr, pr,y, FALSE, NULL);
        }
      }

      if(!isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! x->x) || (! y->x)) 
        {
          if (! (x->x))
            newviewGeneric(tr, pr,x, FALSE, NULL);
          if (! (y->x))
            newviewGeneric(tr, pr,y, FALSE, NULL);
        }
      }				      	

    }

    tr->likelihood = tr->endLH;
  }

  return TRUE;
} 

void restoreTreeFast(tree *tr, partitionList *pr)
{
  removeNodeRestoreBIG(tr, pr, tr->removeNode);
  testInsertRestoreBIG(tr, pr, tr->removeNode, tr->insertNode);
}

static void myfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, stream);

  assert(bytes_written == nmemb);
}

static void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, stream);

  assert(bytes_read == nmemb);
}



/** @brief Write tree to file

    Serialize tree to a file. 

    @todo Document this
*/
static void writeTree(tree *tr, FILE *f)
{
  int 
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    base = tr->nodeBaseAddress;

  myfwrite(&(tr->start->number), sizeof(int), 1, f);
  myfwrite(&base, sizeof(nodeptr), 1, f);
  myfwrite(tr->nodeBaseAddress, sizeof(node), x, f);

}

int ckpCount = 0;

/** @brief Write a checkpoint

    Is checkpoint enabled?

    @todo fill this up
*/
static void writeCheckpoint(tree *tr, partitionList *pr, int state)
{
  int   
    model; 

  char 
    extendedName[2048],
    buf[64];

  FILE 
    *f;

  strcpy(extendedName,  binaryCheckpointName);
  strcat(extendedName, "_");
  sprintf(buf, "%d", ckpCount);
  strcat(extendedName, buf);  

  ckpCount++;

  f = myfopen(extendedName, "w"); 

  /* cdta */   


  tr->ckp.accumulatedTime = accumulatedTime + (gettime() - masterTime);

  tr->ckp.state = state;

  tr->ckp.tr_optimizeRateCategoryInvocations = tr->optimizeRateCategoryInvocations;
  tr->ckp.tr_thoroughInsertion = tr->thoroughInsertion;
  tr->ckp.tr_startLH  = tr->startLH;
  tr->ckp.tr_endLH    = tr->endLH;
  tr->ckp.tr_likelihood = tr->likelihood;
  tr->ckp.tr_bestOfNode = tr->bestOfNode;

  tr->ckp.tr_lhCutoff = tr->lhCutoff;
  tr->ckp.tr_lhAVG    = tr->lhAVG;
  tr->ckp.tr_lhDEC    = tr->lhDEC;     
  tr->ckp.tr_itCount  = tr->itCount;
  tr->ckp.tr_doCutoff = tr->doCutoff;
  /* printf("Acc time: %f\n", tr->ckp.accumulatedTime); */

  /* user stupidity */


  tr->ckp.searchConvergenceCriterion = tr->searchConvergenceCriterion;
  tr->ckp.rateHetModel =  tr->rateHetModel;
  tr->ckp.maxCategories =  tr->maxCategories;
  tr->ckp.NumberOfModels = pr->numberOfPartitions;
  /* tr->ckp.numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1; */
  tr->ckp.originalCrunchedLength = tr->originalCrunchedLength;
  tr->ckp.mxtips = tr->mxtips;
  strcpy(tr->ckp.seq_file, seq_file);

  /* handle user stupidity */


  myfwrite(&(tr->ckp), sizeof(checkPointState), 1, f);

  myfwrite(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfwrite(tr->tree1, sizeof(char), tr->treeStringLength, f);

  myfwrite(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfwrite(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfwrite(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);

  /* need to store this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */

  //TODO: We have to change the way to store the fracchanges
  //myfwrite(tr->fracchanges,  sizeof(double), pr->numberOfPartitions, f);
#if 0 
  myfwrite(&(tr->fracchange),   sizeof(double), 1, f);
#endif


  /* pInfo */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int 
      dataType = pr->partitionData[model]->dataType;

    myfwrite(&(pr->partitionData[model]->numberOfCategories), sizeof(int), 1, f);
    myfwrite(pr->partitionData[model]->perSiteRates, sizeof(double), tr->maxCategories, f);
    myfwrite(pr->partitionData[model]->EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfwrite(pr->partitionData[model]->EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfwrite(pr->partitionData[model]->EI, sizeof(double),  pLengths[dataType].eiLength, f);

    myfwrite(pr->partitionData[model]->frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfwrite(pr->partitionData[model]->tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
    myfwrite(pr->partitionData[model]->substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
    myfwrite(&(pr->partitionData[model]->alpha), sizeof(double), 1, f);
  }



  writeTree(tr, f);

  fclose(f); 

  printBothOpen("\nCheckpoint written to: %s likelihood: %f\n", extendedName, tr->likelihood);
}

static void readTree(tree *tr, partitionList *pr, FILE *f)
{
  int 
    nodeNumber,   
    x = tr->mxtips + 3 * (tr->mxtips - 1);





  nodeptr
    startAddress;

  myfread(&nodeNumber, sizeof(int), 1, f);

  tr->start = tr->nodep[nodeNumber];

  /*printf("Start: %d %d\n", tr->start->number, nodeNumber);*/

  myfread(&startAddress, sizeof(nodeptr), 1, f);

  /*printf("%u %u\n", (size_t)startAddress, (size_t)tr->nodeBaseAddress);*/



  myfread(tr->nodeBaseAddress, sizeof(node), x, f);

  {
    int i;    

    size_t         
      offset;

    boolean 
      addIt;

    if(startAddress > tr->nodeBaseAddress)
    {
      addIt = FALSE;
      offset = (size_t)startAddress - (size_t)tr->nodeBaseAddress;
    }
    else
    {
      addIt = TRUE;
      offset = (size_t)tr->nodeBaseAddress - (size_t)startAddress;
    }       

    for(i = 0; i < x; i++)
    {      	
      if(addIt)
      {	    
        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next + offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back + offset);
      }
      else
      {

        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next - offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back - offset);	   
      } 
    }

  }

  evaluateGeneric(tr, pr, tr->start, TRUE, NULL);

  printBothOpen("RAxML Restart with likelihood: %1.50f\n", tr->likelihood);
}


static void readCheckpoint(tree *tr, partitionList *pr)
{
  int  
    restartErrors = 0,
                  model; 

  FILE 
    *f = myfopen(binaryCheckpointInputName, "r");

  /* cdta */   

  myfread(&(tr->ckp), sizeof(checkPointState), 1, f);



  if(tr->ckp.searchConvergenceCriterion != tr->searchConvergenceCriterion)
  {
    printf("restart error, you are trying to re-start a run where the ML search criterion was turned %s\n", (tr->ckp.searchConvergenceCriterion)?"ON":"OFF");
    restartErrors++;
  }  

  if(tr->ckp.rateHetModel !=  tr->rateHetModel)
  {
    printf("restart error, you are trying to re-start a run with a different model of rate heterogeneity, the checkpoint was obtained under: %s\n", (tr->ckp.rateHetModel == GAMMA)?"GAMMA":"PSR");
    restartErrors++;
  }  

  if(tr->ckp.maxCategories !=  tr->maxCategories)
  {
    printf("restart error, you are trying to re-start a run with %d per-site rate categories, the checkpoint was obtained with: %d\n", tr->maxCategories, tr->ckp.maxCategories);
    restartErrors++;
  }

  if(tr->ckp.NumberOfModels != pr->numberOfPartitions)
  {
    printf("restart error, you are trying to re-start a run with %d partitions, the checkpoint was obtained with: %d partitions\n", (int)pr->numberOfPartitions, tr->ckp.NumberOfModels);
    restartErrors++;      
  }

  /* if(tr->ckp.numBranches != pr->perGeneBranchLengths?pr->numberOfPartitions:1) */
  /* { */
  /*   printf("restart error, you are trying to re-start a run where independent per-site branch length estimates were turned %s\n", (tr->ckp.numBranches > 1)?"ON":"OFF"); */
  /*   restartErrors++; */
  /* } */

  if(tr->ckp.originalCrunchedLength != tr->originalCrunchedLength)
  {
    printf("restart error, you are trying to re-start a run with %d site patterns, the checkpoint was obtained with: %d site patterns\n", tr->ckp.originalCrunchedLength, tr->originalCrunchedLength);
    restartErrors++; 
  }

  if(tr->ckp.mxtips != tr->mxtips)
  {
    printf("restart error, you are trying to re-start a run with %d taxa, the checkpoint was obtained with: %d taxa\n", tr->mxtips, tr->ckp.mxtips);
    restartErrors++; 
  }

  if(strcmp(tr->ckp.seq_file, seq_file) != 0)
  {
    printf("restart error, you are trying to re-start from alignemnt file %s, the checkpoint was obtained with file: %s\n", tr->ckp.seq_file, seq_file);
    restartErrors++; 
  }

  printf("REstart errors: %d\n", restartErrors);

  if(restartErrors > 0)
  {
    printf("User induced errors with the restart from checkpoint, exiting ...\n");

    if(restartErrors > 4)
      printf(" ... maybe you should do field work instead of trying to use a computer ...\n");
    if(restartErrors > 6)
      printf(" ... kala eisai telios ilithios;\n");

    exit(-1);
  }

  tr->ntips = tr->mxtips;

  tr->startLH    = tr->ckp.tr_startLH;
  tr->endLH      = tr->ckp.tr_endLH;
  tr->likelihood = tr->ckp.tr_likelihood;
  tr->bestOfNode = tr->ckp.tr_bestOfNode;

  tr->lhCutoff   = tr->ckp.tr_lhCutoff;
  tr->lhAVG      = tr->ckp.tr_lhAVG;
  tr->lhDEC      = tr->ckp.tr_lhDEC;
  tr->itCount    = tr->ckp.tr_itCount;
  tr->thoroughInsertion       = tr->ckp.tr_thoroughInsertion;



  accumulatedTime = tr->ckp.accumulatedTime;

  /* printf("Accumulated time so far: %f\n", accumulatedTime); */

  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;


  myfread(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfread(tr->tree1, sizeof(char), tr->treeStringLength, f);

  if(tr->searchConvergenceCriterion)
  {
    int bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 0) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 0))
    { 

#ifdef _DEBUG_CHECKPOINTING    
      printf("parsing Tree 0\n");
#endif

      treeReadTopologyString(tr->tree0, tr);   

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }

    bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 1) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 1))
    {

#ifdef _DEBUG_CHECKPOINTING
      printf("parsing Tree 1\n");
#endif

      treeReadTopologyString(tr->tree1, tr); 

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 1, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }
  }

  myfread(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfread(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfread(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);


  /* need to read this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */

  //TODO: Same problem as writing the checkpoint
  //myfread(tr->fracchanges,  sizeof(double), pr->numberOfPartitions, f);
#if 0 
  myfread(&(tr->fracchange),   sizeof(double), 1, f);
#endif

  /* pInfo */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int 
      dataType = pr->partitionData[model]->dataType;

    myfread(&(pr->partitionData[model]->numberOfCategories), sizeof(int), 1, f);
    myfread(pr->partitionData[model]->perSiteRates, sizeof(double), tr->maxCategories, f);
    myfread(pr->partitionData[model]->EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfread(pr->partitionData[model]->EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfread(pr->partitionData[model]->EI, sizeof(double),  pLengths[dataType].eiLength, f);

    myfread(pr->partitionData[model]->frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfread(pr->partitionData[model]->tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
    myfread(pr->partitionData[model]->substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
    myfread(&(pr->partitionData[model]->alpha), sizeof(double), 1, f);

    makeGammaCats(pr->partitionData[model]->alpha, pr->partitionData[model]->gammaRates, 4, tr->useMedian);
  }

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  masterBarrier(THREAD_COPY_INIT_MODEL, tr, pr);
#endif

  updatePerSiteRates(tr, pr, FALSE);

  readTree(tr, pr, f);

  fclose(f); 

}

static void restoreTreeDataValuesFromCheckpoint(tree *tr)
{
  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;  
  tr->thoroughInsertion = tr->ckp.tr_thoroughInsertion;
  tr->likelihood = tr->ckp.tr_likelihood;              
  tr->lhCutoff = tr->ckp.tr_lhCutoff;
  tr->lhAVG    = tr->ckp.tr_lhAVG;
  tr->lhDEC    = tr->ckp.tr_lhDEC;   	 
  tr->itCount = tr->ckp.tr_itCount;
  tr->doCutoff = tr->ckp.tr_doCutoff;
}

void restart(tree *tr, partitionList *pr)
{  
  readCheckpoint(tr, pr);

  switch(tr->ckp.state)
  {
    case REARR_SETTING:      
      break;
    case FAST_SPRS:
      break;
    case SLOW_SPRS:
      break;
    default:
      assert(0);
  }
}

int determineRearrangementSetting(tree *tr, partitionList *pr, analdef *adef, bestlist *bestT, bestlist *bt)
{
  const 
    int MaxFast = 26;

  int 
    i,   
    maxtrav = 5, 
    bestTrav = 5;

  double 
    startLH = tr->likelihood; 

  boolean 
    impr   = TRUE,
           cutoff = tr->doCutoff;

  if(adef->useCheckpoint)
  {
    assert(tr->ckp.state == REARR_SETTING);

    maxtrav = tr->ckp.maxtrav;
    bestTrav = tr->ckp.bestTrav;
    startLH  = tr->ckp.startLH;
    impr     = tr->ckp.impr;      
    cutoff = tr->ckp.cutoff;

    adef->useCheckpoint = FALSE;
  }

  tr->doCutoff = FALSE;      

  resetBestTree(bt);    

#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("MAXTRAV: %d\n", maxtrav);
#endif

  while(impr && maxtrav < MaxFast)
  {	
    recallBestTree(bestT, 1, tr, pr);
    nodeRectifier(tr);                      

    {
      tr->ckp.cutoff = cutoff;	
      tr->ckp.maxtrav = maxtrav;
      tr->ckp.bestTrav = bestTrav;
      tr->ckp.startLH  = startLH;
      tr->ckp.impr = impr;

      writeCheckpoint(tr, pr, REARR_SETTING);
    }

    if (maxtrav > tr->mxtips - 3)  
      maxtrav = tr->mxtips - 3;    

    tr->startLH = tr->endLH = tr->likelihood;

    for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {                	         
      tr->bestOfNode = unlikely;

      if(rearrangeBIG(tr, pr, tr->nodep[i], 1, maxtrav))
      {	     
        if(tr->endLH > tr->startLH)                 	
        {		 	 	      
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;		 
        }	         	       	
      }
    }

    treeEvaluate(tr, pr, 8 ); // 32 * 0.25
    saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("TRAV: %d lh %f\n", maxtrav, tr->likelihood);
#endif

    if(tr->likelihood > startLH)
    {	 
      startLH = tr->likelihood; 	  	  	  
      printLog(tr);	  
      bestTrav = maxtrav;	 
      impr = TRUE;
    }
    else	
      impr = FALSE;	



    if(tr->doCutoff)
    {
      tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));       

      tr->itCount =  tr->itCount + 1;
      tr->lhAVG = 0;
      tr->lhDEC = 0;
    }

    maxtrav += 5;


  }

  recallBestTree(bt, 1, tr, pr);

  tr->doCutoff = cutoff; 

#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("BestTrav %d\n", bestTrav);
#endif

  return bestTrav;     
}





void computeBIGRAPID (tree *tr, partitionList *pr, analdef *adef, boolean estimateModel)
{   
  int
    i,
    impr, 
    bestTrav = 0, 
    rearrangementsMax = 0, 
    rearrangementsMin = 0,    
    thoroughIterations = 0,
    fastIterations = 0;

  double 
    lh = unlikely, 
       previousLh = unlikely, 
       difference, 
       epsilon;              

  bestlist 
    *bestT, 
    *bt;    

  infoList 
    *iList = (infoList*)rax_malloc(sizeof(infoList));

  /* now here is the RAxML hill climbing search algorithm */


  /* initialize two lists of size 1 and size 20 that will keep track of the best 
     and 20 best tree topologies respectively */

  bestT = (bestlist *) rax_malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  bt = (bestlist *) rax_malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  /* initialize an additional data structure used by the search algo, all of this is pretty 
     RAxML-specific and should probably not be in the library */

  initInfoList(iList, 50);

  /* some pretty atbitrary thresholds */

  difference = 10.0;
  epsilon = 0.01;    

  /* Thorough = 0 means that we will do fast SPR inbsertions without optimizing the 
     three branches adjacent to the subtree insertion position via Newton-Raphson 
     */

  tr->thoroughInsertion = FALSE;     

  /* if we are not using a checkpoint and estimateModel is set to TRUE we call the function 
     that optimizes model parameters, such as the CAT model assignment, the alpha paremeter
     or the rates in the GTR matrix. Otherwise we just optimize the branch lengths. Note that 
     the second parameter of treeEvaluate() controls how many times we will iterate over all branches 
     of the tree until we give up, provided that, the br-len opt. has not converged before.
     */

  if(!adef->useCheckpoint)
  {
    if(estimateModel)
      modOpt(tr, pr, 10.0);
    else
      treeEvaluate(tr, pr, 64); // 32 * 2
  }

  /* print some stuff to the RAxML_log file */

  printLog(tr); 

  /* save the current tree (which is the input tree parsed via -t in the bestT list */

  saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  /* if the rearrangmenet radius has been set by the user ie. adef->initailSet == TRUE 
     then just set the apppropriate parameter.
     Otherwise, call the function  determineRearrangementSetting() that seeks 
     for the best radius by executing SPR moves on the initial tree with different radii
     and returns the smallest radius that yields the best log likelihood score after 
     applying one cycle of SPR moves to the tree 
     */

  if(!adef->initialSet)   
  {
    if((!adef->useCheckpoint) || (adef->useCheckpoint && tr->ckp.state == REARR_SETTING))
    {
      bestTrav = adef->bestTrav = determineRearrangementSetting(tr, pr, adef, bestT, bt);
      printBothOpen("\nBest rearrangement radius: %d\n", bestTrav);
    }
  }
  else
  {
    bestTrav = adef->bestTrav = adef->initial;       
    printBothOpen("\nUser-defined rearrangement radius: %d\n", bestTrav);
  }


  /* some checkpointing noise */
  if(!(adef->useCheckpoint && (tr->ckp.state == FAST_SPRS || tr->ckp.state == SLOW_SPRS)))
  {      

    /* optimize model params more thoroughly or just optimize branch lengths */
    if(estimateModel)
      modOpt(tr, pr, 5.0);
    else
      treeEvaluate(tr, pr, 32);   // 32 * 1
  }

  /* save the current tree again, while the topology has not changed, the branch lengths have changed in the meantime, hence
     we need to store them again */

  saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  /* set the loop variable to TRUE */

  impr = 1;

  /* this is for the additional RAxML heuristics described imn this paper here:

     A. Stamatakis,  F. Blagojevic, C.D. Antonopoulos, D.S. Nikolopoulos: "Exploring new Search Algorithms and Hardware for Phylogenetics: RAxML meets the IBM Cell". 
     In Journal of VLSI Signal Processing Systems, 48(3):271-286, 2007.

     This is turned on by default 
     */


  if(tr->doCutoff)
    tr->itCount = 0;

  /* figure out where to continue computations if we restarted from a checkpoint */

  if(adef->useCheckpoint && tr->ckp.state == FAST_SPRS)
    goto START_FAST_SPRS;

  if(adef->useCheckpoint && tr->ckp.state == SLOW_SPRS)
    goto START_SLOW_SPRS;

  while(impr)
  {
START_FAST_SPRS:
    /* if re-starting from checkpoint set the required variable values to the 
       values that they had when the checkpoint was written */

    if(adef->useCheckpoint && tr->ckp.state == FAST_SPRS)
    {	    
      impr = tr->ckp.impr;	  
      bestTrav = tr->ckp.bestTrav;	  
      rearrangementsMax = tr->ckp.rearrangementsMax;
      rearrangementsMin = tr->ckp.rearrangementsMin;
      thoroughIterations = tr->ckp.thoroughIterations;
      fastIterations = tr->ckp.fastIterations;  
      lh = tr->ckp.lh;
      previousLh = tr->ckp.previousLh;
      difference = tr->ckp.difference;
      epsilon    = tr->ckp.epsilon;  

      restoreTreeDataValuesFromCheckpoint(tr);	  

      adef->useCheckpoint = FALSE;
    }
    else
      /* otherwise, restore the currently best tree */
      recallBestTree(bestT, 1, tr, pr);

    /* save states of algorithmic/heuristic variables for printing the next checkpoint */


    tr->ckp.impr = impr;	
    tr->ckp.bestTrav = bestTrav;      
    tr->ckp.rearrangementsMax = rearrangementsMax;
    tr->ckp.rearrangementsMin = rearrangementsMin;
    tr->ckp.thoroughIterations = thoroughIterations;
    tr->ckp.fastIterations = fastIterations;  
    tr->ckp.lh = lh;
    tr->ckp.previousLh = previousLh;
    tr->ckp.difference = difference;
    tr->ckp.epsilon    = epsilon;              
    tr->ckp.bestTrav = bestTrav;       
    tr->ckp.impr = impr;                 

    /* write a binary checkpoint */
    writeCheckpoint(tr, pr, FAST_SPRS);

    /* this is the aforementioned convergence criterion that requires computing the RF,
       let's not worry about this right now */

    if(tr->searchConvergenceCriterion)
    {
      int bCounter = 0;	  	      	 	  	  	

      if(fastIterations > 1)
        cleanupHashTable(tr->h, (fastIterations % 2));		
      
      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, fastIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE, tr->threadID);	    

      {
        char 
          *buffer = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
#ifdef _DEBUG_CHECKPOINTING
        printf("Storing tree in slot %d\n", fastIterations % 2);
#endif

        Tree2String(buffer, tr, pr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);

        if(fastIterations % 2 == 0)	      
          memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
        else
          memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    

        rax_free(buffer);
      }


      assert(bCounter == tr->mxtips - 3);	    	   

      if(fastIterations > 0)
      {
        double rrf = convergenceCriterion(tr->h, tr->mxtips);

        if(rrf <= 0.01) /* 1% cutoff */
        {
          printBothOpen("ML fast search converged at fast SPR cycle %d with stopping criterion\n", fastIterations);
          printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
          cleanupHashTable(tr->h, 0);
          cleanupHashTable(tr->h, 1);
          goto cleanup_fast;
        }
        else		    
          printBothOpen("ML search convergence criterion fast cycle %d->%d Relative Robinson-Foulds %f\n", fastIterations - 1, fastIterations, rrf);
      }
    }


    /* count how many fast iterations with so-called fast SPR moves we have executed */

    fastIterations++;	

    /* optimize branch lengths */

    treeEvaluate(tr, pr, 32);  // 32 * 1 = 32

    /* save the tree with those branch lengths again */

    saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

    /* print the log likelihood */

    printLog(tr);    

    /* print this intermediate tree to file */

    printResult(tr, pr, adef, FALSE);

    /* update the current best likelihood */

    lh = previousLh = tr->likelihood;

    /* in here we actually do a cycle of SPR moves */

    treeOptimizeRapid(tr, pr, 1, bestTrav, adef, bt, iList);

    /* set impr to 0 since in the immediately following for loop we check if the SPR moves above have generated 
       a better tree */

    impr = 0;

    /* loop over the 20 best trees generated by the fast SPR moves, and check if they improve the likelihood after all of their branch lengths
       have been optimized */

    for(i = 1; i <= bt->nvalid; i++)
    {	    	
      /* restore tree i from list generated by treeOptimizeRapid */

      recallBestTree(bt, i, tr, pr);

      /* optimize branch lengths of this tree */

      treeEvaluate(tr, pr, 8); // 0.25 * 32

      /* calc. the likelihood improvement */

      difference = ((tr->likelihood > previousLh)? 
          tr->likelihood - previousLh: 
          previousLh - tr->likelihood); 	    

      /* if the likelihood has improved save the current tree as best tree and continue */
      /* note that we always compre this tree to the likelihood of the previous best tree */

      if(tr->likelihood > lh && difference > epsilon)
      {
        impr = 1;	       
        lh = tr->likelihood;	       	     
        saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

      }	   	   
    }
#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("FAST LH: %f\n", lh);
#endif


  }

  /* needed for this RF-based convergence criterion that I actually describe in here:

     A. Stamatakis: "Phylogenetic Search Algorithms for Maximum Likelihood". In M. Elloumi, A.Y. Zomaya, editors. 
     Algorithms in Computational Biology: techniques, Approaches and Applications, John Wiley and Sons

     a copy of this book is in my office */

  if(tr->searchConvergenceCriterion)
  {
    cleanupHashTable(tr->h, 0);
    cleanupHashTable(tr->h, 1);
  }

cleanup_fast:  
  /*
     now we have jumped out of the loop that executes 
     fast SPRs, and next we will execute a loop that executes throough SPR cycles (with SPR moves 
     that optimize via newton-Raphson all adjacent branches to the insertion point) 
     until no through SPR move can be found that improves the likelihood further. A classic 
     hill climbing algo.
     */

  tr->thoroughInsertion = TRUE;
  impr = 1;

  /* restore the currently best tree. this si actually required, because we do not know which tree
     is actually stored in the tree data structure when the above loop exits */

  recallBestTree(bestT, 1, tr, pr);

  {
    /* RE-TRAVERSE THE ENTIRE TREE */

    evaluateGeneric(tr, pr, tr->start, TRUE, NULL);
#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("After Fast SPRs Final %f\n", tr->likelihood);   
#endif
  }

  /* optimize model params (including branch lengths) or just 
     optimize branch lengths and leave the other model parameters (GTR rates, alhpa) 
     alone */

  if(estimateModel)
    modOpt(tr, pr, 1.0);
  else
    treeEvaluate(tr, pr, 32 ); //32 * 1

  /* start loop that executes thorough SPR cycles */

  while(1)
  {	 
    /* once again if we want to restart from a checkpoint that was written during this loop we need
       to restore the values of the variables appropriately */
START_SLOW_SPRS:
    if(adef->useCheckpoint && tr->ckp.state == SLOW_SPRS)
    {	        
      impr = tr->ckp.impr;	 
      bestTrav = tr->ckp.bestTrav;	 
      rearrangementsMax = tr->ckp.rearrangementsMax;
      rearrangementsMin = tr->ckp.rearrangementsMin;
      thoroughIterations = tr->ckp.thoroughIterations;
      fastIterations = tr->ckp.fastIterations;     
      lh = tr->ckp.lh;
      previousLh = tr->ckp.previousLh;
      difference = tr->ckp.difference;
      epsilon    = tr->ckp.epsilon;                    

      restoreTreeDataValuesFromCheckpoint(tr);	  

      adef->useCheckpoint = FALSE;
    }
    else
      /* otherwise we restore the currently best tree and load it from bestT into our tree data 
         structuire tr */
      recallBestTree(bestT, 1, tr, pr);

    /* now, we write a checkpoint */

    tr->ckp.impr = impr;      
    tr->ckp.bestTrav = bestTrav;      
    tr->ckp.rearrangementsMax = rearrangementsMax;
    tr->ckp.rearrangementsMin = rearrangementsMin;
    tr->ckp.thoroughIterations = thoroughIterations;
    tr->ckp.fastIterations = fastIterations;	  
    tr->ckp.lh = lh;
    tr->ckp.previousLh = previousLh;
    tr->ckp.difference = difference;
    tr->ckp.epsilon    = epsilon;              
    tr->ckp.bestTrav = bestTrav;       
    tr->ckp.impr = impr;                 

    /* write binary checkpoint to file */
    writeCheckpoint(tr, pr, SLOW_SPRS);

    if(impr)
    {
      /* if the logl has improved write out some stuff and adapt the rearrangement radii */
      printResult(tr, pr, adef, FALSE);
      /* minimum rearrangement radius */
      rearrangementsMin = 1;
      /* max radius, this is probably something I need to explain at the whiteboard */
      rearrangementsMax = adef->stepwidth;	

      /* once again the convergence criterion */

      if(tr->searchConvergenceCriterion)
      {
        int bCounter = 0;	      

        if(thoroughIterations > 1)
          cleanupHashTable(tr->h, (thoroughIterations % 2));		

        bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, thoroughIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
            &bCounter, 1, FALSE, FALSE, tr->threadID);	    


        {
          char 
            *buffer = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));

#ifdef _DEBUG_CHECKPOINTING		
          printf("Storing tree in slot %d\n", thoroughIterations % 2);
#endif

          Tree2String(buffer, tr, pr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);

          if(thoroughIterations % 2 == 0)	      
            memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
          else
            memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    

          rax_free(buffer);
        }

        assert(bCounter == tr->mxtips - 3);

        if(thoroughIterations > 0)
        {
          double rrf = convergenceCriterion(tr->h, tr->mxtips);

          if(rrf <= 0.01) /* 1% cutoff */
          {
            printBothOpen("ML search converged at thorough SPR cycle %d with stopping criterion\n", thoroughIterations);
            printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
            goto cleanup;
          }
          else		    
            printBothOpen("ML search convergence criterion thorough cycle %d->%d Relative Robinson-Foulds %f\n", thoroughIterations - 1, thoroughIterations, rrf);
        }
      }



      thoroughIterations++;	  
    }			  			
    else
    {

      /* if the lnl has not imrpved by the current SPR cycle adapt the min and max rearrangemnt radii and try again */

      rearrangementsMax += adef->stepwidth;
      rearrangementsMin += adef->stepwidth; 	        	      

      /* if we have already tried them then abandon this loop, the search has converged */
      if(rearrangementsMax > adef->max_rearrange)	     	     	 
        goto cleanup; 	   
    }

    /* optimize branch lengths of best tree */

    treeEvaluate(tr, pr, 32 ); // 32 * 1

    /* do some bokkeeping and printouts again */
    previousLh = lh = tr->likelihood;	      
    saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
    printLog(tr);

    /* do a cycle of thorough SPR moves with the minimum and maximum rearrangement radii */

    treeOptimizeRapid(tr, pr, rearrangementsMin, rearrangementsMax, adef, bt, iList);

    impr = 0;			      		            

    /* once again get the best 20 trees produced by the SPR cycle, load them from the bt tree list into tr
       optimize their branch lengths and figure out if the LnL of the tree has improved */

    for(i = 1; i <= bt->nvalid; i++)
    {		 
      recallBestTree(bt, i, tr, pr);

      treeEvaluate(tr, pr, 8); // 0.25	* 32

      difference = ((tr->likelihood > previousLh)? 
          tr->likelihood - previousLh: 
          previousLh - tr->likelihood); 	    
      if(tr->likelihood > lh && difference > epsilon)
      {
        impr = 1;	       
        lh = tr->likelihood;	  	     
        saveBestTree(bestT, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
      }	   	   
    }  

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("SLOW LH: %f\n", lh);              
#endif
  }

cleanup: 

  /* do a final full tree traversal, not sure if this is required here */

  {
    evaluateGeneric(tr, pr, tr->start, TRUE, NULL);

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("After SLOW SPRs Final %f\n", tr->likelihood);   
#endif
  }

  /* free data structures */

  if(tr->searchConvergenceCriterion)
  {
    freeBitVectors(tr->bitVectors, 2 * tr->mxtips);
    rax_free(tr->bitVectors);
    freeHashTable(tr->h);
    rax_free(tr->h);
  }

  freeBestTree(bestT);
  rax_free(bestT);
  freeBestTree(bt);
  rax_free(bt);
  freeInfoList(iList);  
  rax_free(iList);

  printLog(tr);

  printResult(tr, pr, adef, TRUE);

  /* and we are done, return to main() in axml.c  */

}



/* The number of maximum smoothing iterations is given explicitely */
boolean 
treeEvaluate (tree *tr, partitionList *pr, int maxSmoothIterations)       /* Evaluate a user tree */
{
  smoothTree(tr, pr, maxSmoothIterations); /* former (32 * smoothFactor) */

  evaluateGeneric(tr, pr, tr->start, FALSE, NULL);

  return TRUE;
}

/** @brief Perform an NNI move

    Modify the topology of tree \a tr by performing an NNI (Neighbour Neighbor
    Interchange) move at node \a p. Perform one of the two possible NNI moves
    based on whether \a swap is set to 1 or 2.
*/
void NNI(tree * tr, nodeptr p, int swap)
{
  nodeptr       q, tmp;

  q = p->back;
  assert(!isTip(q->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));


  if(swap == 1)
   {
     tmp = p->next->back;
     hookupFull(p->next, q->next->back, q->next->z);
     hookupFull(q->next, tmp,           p->next->z);
   }
  else
   {
      tmp = p->next->next->back;
      hookupFull(p->next->next, q->next->back, q->next->z);
      hookupFull(q->next,       tmp,           p->next->next->z);
   }
}
