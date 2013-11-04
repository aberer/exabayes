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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include "mem_alloc.h"

#ifndef WIN32
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

static const double MNBRAK_GOLD =    1.618034;
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =      1.e-5;
static const double BRENT_CGOLD =   0.3819660;

extern int optimizeRatesInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[20];




/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


/* the following function is used to set rates in the Q matrix 
   the data structure called symmetryVector is used to 
   define the symmetries between rates as they are specified 
   in some of the secondary structure substitution models that 
   generally don't use GTR matrices but more restricted forms thereof */

static void setRateModel(partitionList *pr, int model, double rate, int position)
{
  int
    states   = pr->partitionData[model]->states,
    numRates = (states * states - states) / 2;

  if(pr->partitionData[model]->dataType == DNA_DATA)
    assert(position >= 0 && position < (numRates - 1));
  else
    assert(position >= 0 && position < numRates);

  assert(pr->partitionData[model]->dataType != BINARY_DATA);

  assert(rate >= RATE_MIN && rate <= RATE_MAX);

  if(pr->partitionData[model]->nonGTR)
    {    
      int 
	i, 
	k = pr->partitionData[model]->symmetryVector[position];

      assert(pr->partitionData[model]->dataType == SECONDARY_DATA ||
	     pr->partitionData[model]->dataType == SECONDARY_DATA_6 ||
	     pr->partitionData[model]->dataType == SECONDARY_DATA_7);

      if(k == -1)
	pr->partitionData[model]->substRates[position] = 0.0;
      else
	{
	  if(k == pr->partitionData[model]->symmetryVector[numRates - 1])
	    {
	      for(i = 0; i < numRates - 1; i++)
		if(pr->partitionData[model]->symmetryVector[i] == k)
		  pr->partitionData[model]->substRates[position] = 1.0;
	    }
	  else
	    {
	      for(i = 0; i < numRates - 1; i++)
		{
		  if(pr->partitionData[model]->symmetryVector[i] == k)
		    pr->partitionData[model]->substRates[i] = rate;
		}	      	     
	    }
	}
    }
  else
    pr->partitionData[model]->substRates[position] = rate;
}


/* 
   the following three functions are used to link/unlink parameters 
   between partitions. This should work in a generic way, however 
   this is so far mainly used for linking unlinking GTR matrix parameter 
   estimates across different protein data partitions.
   Generally this mechanism can also be used for linking/inlinking alpha paremeters 
   between partitions and the like.
   However, all alpha parameter estimates for all partitions and GTR estimates for 
   DNA partitions are unlinked by default. This is actually hard-coded 
   in here. 
*/

/* initializwe a parameter linkage list for a certain parameter type (can be whatever).
   the input is an integer vector that contaions NumberOfModels (numberOfPartitions) elements.

   if we want to have all alpha parameters unlinked and have say 4 partitions the input 
   vector would look like this: {0, 1, 2, 3}, if we want to link partitions 0 and 3 the vector 
   should look like this: {0, 1, 2, 0} 
*/

static linkageList* initLinkageList(int *linkList, partitionList *pr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  
  linkageList 
    *ll = (linkageList*)rax_malloc(sizeof(linkageList));
    
  /* figure out how many distinct parameters we need to estimate 
     in total, if all parameters are linked the result will be 1 if all 
     are unlinked the result will be pr->numberOfPartitions */
  
  for(i = 0; i < pr->numberOfPartitions; i++)
    {
      if(linkList[i] > numberOfModels)
	numberOfModels = linkList[i];
    }

  numberOfModels++;
  
  /* allocate the linkage list data structure that containes information which parameters of which partition are 
     linked with each other.

     Note that we need a separate invocation of initLinkageList() and a separate linkage list 
     for each parameter type */

  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)rax_malloc(sizeof(linkageData) * numberOfModels);

  /* noe loop over the number of free parameters and assign the corresponding partitions to each parameter */

  for(i = 0; i < numberOfModels; i++)
    {
      /* 
	 the valid flag is used for distinguishing between DNA and protein data partitions.
	 This can be used to enable/disable parameter optimization for the paremeter 
	 associated to the corresponding partitions. This deature is used in optRatesGeneric 
	 to first optimize all DNA GTR rate matrices and then all PROT GTR rate matrices */

      ll->ld[i].valid = TRUE;
      partitions = 0;

      /* now figure out how many partitions share this joint parameter */

      for(k = 0; k < pr->numberOfPartitions; k++)
	if(linkList[k] == i)
	  partitions++;	    

      /* assign a list to store the partitions that share the parameter */

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)rax_malloc(sizeof(int) * partitions);
      
      /* now store the respective partition indices in this list */
      
      for(k = 0, pos = 0; k < pr->numberOfPartitions; k++)
	if(linkList[k] == i)
	  ll->ld[i].partitionList[pos++] = k;
    }

  /* return the linkage list for the parameter */

  return ll;
}

/* dedicated helper function to initialize the linkage list, that is, essentiaylly compute 
   the integer vector int *linkList used above for linking GTR models.
   
   Once again, this is hard-coded in RAxML, because users can not influence the linking.

*/
   

static linkageList* initLinkageListGTR(partitionList *pr)
{
  int
    i,
    *links = (int*)rax_malloc(sizeof(int) * pr->numberOfPartitions),
    firstAA = pr->numberOfPartitions + 2,
    countGTR = 0,
    countOtherModel = 0;
  
  linkageList
    * ll;

  /* here we only want to figure out if either all prot data partitions 
     are supposed to use a joint GTR prot subst matrix or not 

    We either allow ALL prot partitions to use a shared/joint estimate of the GTR matrix or not,
    things like having one prot partition evolving under WAG and the others under a joint GTR estimate are 
    not allowed.
  */

  for(i = 0; i < pr->numberOfPartitions; i++)
    {     
      if(pr->partitionData[i]->dataType == AA_DATA)
	{
	  if(pr->partitionData[i]->protModels == GTR)
	    {
	      if(i < firstAA)
		firstAA = i;
	      countGTR++;
	    }
	  else
	    countOtherModel++;
	}
    }
  
  assert((countGTR > 0 && countOtherModel == 0) || (countGTR == 0 && countOtherModel > 0) ||  (countGTR == 0 && countOtherModel == 0));

  /* if there is no joint GTR matrix optimization for protein data partitions we can unlink rate matrix calculations for all partitions */

  if(countGTR == 0)
    {
      for(i = 0; i < pr->numberOfPartitions; i++)
	links[i] = i;
    }
  else
    {
      /* otherwise we let all partitions, except for the protein partitions use 
	 unlinked rate matrices while we link the GTR rate matrices of all 
	 protein data partitions */
      for(i = 0; i < pr->numberOfPartitions; i++)
	{
	  switch(pr->partitionData[i]->dataType)
	    {	   
	    case DNA_DATA:
	    case BINARY_DATA:
	    case GENERIC_32:
	    case GENERIC_64:
	    case SECONDARY_DATA:
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7: 
	      links[i] = i;
	      break;
	    case AA_DATA:	  
	      links[i] = firstAA;
	      break;
	    default:
	      assert(0);
	    }
	}
    }
  

  /* we can now pass an appropriate integer vector to the linkage list initialization function :-) */

  ll = initLinkageList(links, pr);

  rax_free(links);
  
  return ll;
}

/* free linkage list data structure */

static void freeLinkageList( linkageList* ll)
{
  int i;    

  for(i = 0; i < ll->entries; i++)    
    rax_free(ll->ld[i].partitionList);         

  rax_free(ll->ld);
  rax_free(ll);   
}

#define ALPHA_F 0
#define INVAR_F 1
#define RATE_F  2


/* function that evaluates the change to a parameter */

static void evaluateChange(tree *tr, partitionList *pr, int rateNumber, double *value, double *result, boolean* converged, int whichFunction, int numberOfModels, linkageList *ll)
{ 
  int i, k, pos;

  switch(whichFunction)
    {    
    case RATE_F:
      /* loop over linkage list entries for Q matrix rates */

      for(i = 0, pos = 0; i < ll->entries; i++)
	{
	  /* if valid, i.e. we want to assess the change for a certain partition data type, e.g., DNA */
	  if(ll->ld[i].valid)
	    {
	      /* if the iterative numerical procedure for the partitions sharing the parameter 
		 has already converged don't do anything.
		 This stuff is required for parallel load balanicing as described in:
		 
		 A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009, Vienna, Austria, September 2009.
	      */
	      if(converged[pos])
		{		 
		  for(k = 0; k < ll->ld[i].partitions; k++)
		    pr->partitionData[ll->ld[i].partitionList[k]]->executeModel = FALSE;
		}
	      else
		{
		  /* loop over partitions that share the same parameter */

		  for(k = 0; k < ll->ld[i].partitions; k++)
		    {
		      int index = ll->ld[i].partitionList[k];		  
		      
		      /* update them to new, proposed value for which we want to obtain the likelihood */
	      
		      setRateModel(pr, index, value[pos], rateNumber);

		      /* in the case of rates, i.e., Q matrix, re-compute the corresponding eigenvector eigenvalue decomposition */

		      initReversibleGTR(tr, pr, index);
		    }
		}
	      pos++;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
	    	  pr->partitionData[ll->ld[i].partitionList[k]]->executeModel = FALSE;
	    }
	 
	}

      assert(pos == numberOfModels);

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))      
      masterBarrier(THREAD_OPT_RATE, tr, pr);
#else
      /* and compute the likelihood by doing a full tree traversal :-) */
      evaluateGeneric(tr, pr, tr->start, TRUE);
#endif     
      
      /* update likelihoods and the sum of per-partition likelihoods for those partitions that share the parameter.
	 that's the likelihood we actually want to optimize here */

      for(i = 0, pos = 0; i < ll->entries; i++)	
	{
	  if(ll->ld[i].valid)
	    {
	      result[pos] = 0.0;
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int index = ll->ld[i].partitionList[k];

		  assert(pr->partitionData[index]->partitionLH <= 0.0);

		  result[pos] -= pr->partitionData[index]->partitionLH;
		  
		}
	      pos++;
	    }
	   for(k = 0; k < ll->ld[i].partitions; k++)
	     {
	       int index = ll->ld[i].partitionList[k];
	       pr->partitionData[index]->executeModel = TRUE;
	     }	  
	}

      assert(pos == numberOfModels);
      break;
    case ALPHA_F:

      /* analoguos to the rate stuff above */

      for(i = 0; i < ll->entries; i++)
	{
	  if(converged[i])
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
	    	  pr->partitionData[ll->ld[i].partitionList[k]]->executeModel = FALSE;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int index = ll->ld[i].partitionList[k];
		  pr->partitionData[index]->executeModel = TRUE;
		  pr->partitionData[index]->alpha = value[i];

		  /* re-compute the discrete gamma function approximation for the new alpha parameter */
		  makeGammaCats(pr->partitionData[index]->alpha, pr->partitionData[index]->gammaRates, 4, tr->useMedian);
		}
	    }
	}
#if (defined( _USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
      masterBarrier(THREAD_OPT_ALPHA, tr, pr);
#else  
      evaluateGeneric(tr, pr, tr->start, TRUE);
#endif
            
      for(i = 0; i < ll->entries; i++)	
	{	  
	  result[i] = 0.0;
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {
	      int index = ll->ld[i].partitionList[k];
	      	      
	      assert(pr->partitionData[index]->partitionLH <= 0.0);
	      
	      result[i] -= pr->partitionData[index]->partitionLH;
	      pr->partitionData[index]->executeModel = TRUE;
	    }
	}
      break;
    default:
      assert(0);	
    }

}

/* generic implementation of Brent's algorithm for one-dimensional parameter optimization */

static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, tree *tr, partitionList *pr, linkageList *ll, double lim_inf, double lim_sup)
{
  int iter, i;
  double 
    *a     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *b     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *d     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *etemp = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fu    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fv    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fw    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fx    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *p     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *q     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *r     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *tol1  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *tol2  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *u     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *v     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *w     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *x     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *xm    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *e     = (double *)rax_malloc(sizeof(double) * numberOfModels);
  boolean *converged = (boolean *)rax_malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;
  
  for(i = 0; i < numberOfModels; i++)    
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      e[i] = 0.0;
      d[i] = 0.0;
    }

  for(i = 0; i < numberOfModels; i++)
    {
      a[i]=((ax[i] < cx[i]) ? ax[i] : cx[i]);
      b[i]=((ax[i] > cx[i]) ? ax[i] : cx[i]);
      x[i] = w[i] = v[i] = bx[i];
      fw[i] = fv[i] = fx[i] = fb[i];
    }

  for(i = 0; i < numberOfModels; i++)
    {      
      assert(a[i] >= lim_inf && a[i] <= lim_sup);
      assert(b[i] >= lim_inf && b[i] <= lim_sup);
      assert(x[i] >= lim_inf && x[i] <= lim_sup);
      assert(v[i] >= lim_inf && v[i] <= lim_sup);
      assert(w[i] >= lim_inf && w[i] <= lim_sup);
    }
  
  

  for(iter = 1; iter <= ITMAX; iter++)
    {
      allConverged = TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
	allConverged = allConverged && converged[i];

      if(allConverged)
	{
	  rax_free(converged);
	  rax_free(a);
	  rax_free(b);
	  rax_free(d);
	  rax_free(etemp);
	  rax_free(fu);
	  rax_free(fv);
	  rax_free(fw);
	  rax_free(fx);
	  rax_free(p);
	  rax_free(q);
	  rax_free(r);
	  rax_free(tol1);
	  rax_free(tol2);
	  rax_free(u);
	  rax_free(v);
	  rax_free(w);
	  rax_free(x);
	  rax_free(xm);
	  rax_free(e);
	  return;
	}     

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {	     	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
  
	      xm[i] = 0.5 * (a[i] + b[i]);
	      tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BRENT_ZEPS);
	  
	      if(fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
		{		 
		  result[i] =  -fx[i];
		  xmin[i]   = x[i];
		  converged[i] = TRUE;		  
		}
	      else
		{
		  if(fabs(e[i]) > tol1[i])
		    {		     
		      r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
		      q[i] = (x[i] - v[i]) * (fx[i] - fw[i]);
		      p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
		      q[i] = 2.0 * (q[i] - r[i]);
		      if(q[i] > 0.0)
			p[i] = -p[i];
		      q[i] = fabs(q[i]);
		      etemp[i] = e[i];
		      e[i] = d[i];
		      if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i])) || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
			d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
		      else
			{
			  d[i] = p[i] / q[i];
			  u[i] = x[i] + d[i];
			  if( u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
			    d[i] = SIGN(tol1[i], xm[i] - x[i]);
			}
		    }
		  else
		    {		     
		      d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i]: b[i] - x[i]));
		    }
		  u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]): (x[i] +SIGN(tol1[i], d[i])));
		}

	      if(!converged[i])
		assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
                 
      evaluateChange(tr, pr, rateNumber, u, fu, converged, whichFunction, numberOfModels, ll);

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {
	      if(fu[i] <= fx[i])
		{
		  if(u[i] >= x[i])
		    a[i] = x[i];
		  else
		    b[i] = x[i];
		  
		  SHFT(v[i],w[i],x[i],u[i]);
		  SHFT(fv[i],fw[i],fx[i],fu[i]);
		}
	      else
		{
		  if(u[i] < x[i])
		    a[i] = u[i];
		  else
		    b[i] = u[i];
		  
		  if(fu[i] <= fw[i] || w[i] == x[i])
		    {
		      v[i] = w[i];
		      w[i] = u[i];
		      fv[i] = fw[i];
		      fw[i] = fu[i];
		    }
		  else
		    {
		      if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
			{
			  v[i] = u[i];
			  fv[i] = fu[i];
			}
		    }	    
		}
	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
	      assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
    }

  rax_free(converged);
  rax_free(a);
  rax_free(b);
  rax_free(d);
  rax_free(etemp);
  rax_free(fu);
  rax_free(fv);
  rax_free(fw);
  rax_free(fx);
  rax_free(p);
  rax_free(q);
  rax_free(r);
  rax_free(tol1);
  rax_free(tol2);
  rax_free(u);
  rax_free(v);
  rax_free(w);
  rax_free(x);
  rax_free(xm);
  rax_free(e);

  printf("\n. Too many iterations in BRENT !");
  assert(0);
}

/* generic bracketing function required for Brent's algorithm. For details please see the corresponding chapter in the book Numerical Recipees in C */

static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double lim_inf, double lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, partitionList *pr, linkageList *ll)
{
  double 
    *ulim = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *u    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *r    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *q    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fu   = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *dum  = (double *)rax_malloc(sizeof(double) * numberOfModels), 
    *temp = (double *)rax_malloc(sizeof(double) * numberOfModels);
  
  int 
    i,
    *state    = (int *)rax_malloc(sizeof(int) * numberOfModels),
    *endState = (int *)rax_malloc(sizeof(int) * numberOfModels);

  boolean *converged = (boolean *)rax_malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      state[i] = 0;
      endState[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > lim_sup) 	
	param[i] = ax[i] = lim_sup;
      
      if(param[i] < lim_inf) 
	param[i] = ax[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
   
  
  evaluateChange(tr, pr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup) 
	param[i] = bx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = bx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
  evaluateChange(tr, pr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
	{	  
	  SHFT(dum[i],ax[i],bx[i],dum[i]);
	  SHFT(dum[i],fa[i],fb[i],dum[i]);
	}
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup) 
	param[i] = cx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = cx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
 
  evaluateChange(tr, pr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll);

   while(1) 
     {       
       allConverged = TRUE;

       for(i = 0; i < numberOfModels && allConverged; i++)
	 allConverged = allConverged && converged[i];

       if(allConverged)
	 {
	   for(i = 0; i < numberOfModels; i++)
	     {	       
	       if(ax[i] > lim_sup) 
		 ax[i] = lim_sup;
	       if(ax[i] < lim_inf) 
		 ax[i] = lim_inf;

	       if(bx[i] > lim_sup) 
		 bx[i] = lim_sup;
	       if(bx[i] < lim_inf) 
		 bx[i] = lim_inf;
	       
	       if(cx[i] > lim_sup) 
		 cx[i] = lim_sup;
	       if(cx[i] < lim_inf) 
		 cx[i] = lim_inf;
	     }

	   rax_free(converged);
	   rax_free(ulim);
	   rax_free(u);
	   rax_free(r);
	   rax_free(q);
	   rax_free(fu);
	   rax_free(dum); 
	   rax_free(temp);
	   rax_free(state);   
	   rax_free(endState);
	   return 0;
	   
	 }

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {
	       switch(state[i])
		 {
		 case 0:
		   endState[i] = 0;
		   if(!(fb[i] > fc[i]))		         
		     converged[i] = TRUE;		       		     
		   else
		     {
		   
		       if(ax[i] > lim_sup) 
			 ax[i] = lim_sup;
		       if(ax[i] < lim_inf) 
			 ax[i] = lim_inf;
		       if(bx[i] > lim_sup) 
			 bx[i] = lim_sup;
		       if(bx[i] < lim_inf) 
			 bx[i] = lim_inf;
		       if(cx[i] > lim_sup) 
			 cx[i] = lim_sup;
		       if(cx[i] < lim_inf) 
			 cx[i] = lim_inf;
		       
		       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
		       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
		       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
			 (2.0*SIGN(MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
		       
		       ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
		       
		       if(u[i] > lim_sup) 
			 u[i] = lim_sup;
		       if(u[i] < lim_inf) 
			 u[i] = lim_inf;
		       if(ulim[i] > lim_sup) 
			 ulim[i] = lim_sup;
		       if(ulim[i] < lim_inf) 
			 ulim[i] = lim_inf;
		       
		       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
			 {
			   param[i] = u[i];
			   if(param[i] > lim_sup) 			     
			     param[i] = u[i] = lim_sup;
			   if(param[i] < lim_inf)
			     param[i] = u[i] = lim_inf;
			   endState[i] = 1;
			 }
		       else 
			 {
			   if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
			     {
			       param[i] = u[i];
			       if(param[i] > lim_sup) 
				 param[i] = u[i] = lim_sup;
			       if(param[i] < lim_inf) 
				 param[i] = u[i] = lim_inf;
			       endState[i] = 2;
			     }		  	       
			   else
			     {
			       if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
				 {
				   u[i] = ulim[i];
				   param[i] = u[i];	
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = ulim[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = ulim[i] = lim_inf;
				   endState[i] = 0;
				 }		  		
			       else 
				 {		  
				   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
				   param[i] = u[i];
				   endState[i] = 0;
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = lim_inf;
				 }
			     }	  
			 }
		     }
		   break;
		 case 1:
		   endState[i] = 0;
		   break;
		 case 2:
		   endState[i] = 3;
		   break;
		 default:
		   assert(0);
		 }
	       assert(param[i] >= lim_inf && param[i] <= lim_sup);
	     }
	 }
             
       evaluateChange(tr, pr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll);

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {	       
	       switch(endState[i])
		 {
		 case 0:
		   fu[i] = temp[i];
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 case 1:
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {
		       ax[i]=(bx[i]);
		       bx[i]=u[i];
		       fa[i]=(fb[i]);
		       fb[i]=fu[i]; 
		       converged[i] = TRUE;		      
		     } 
		   else 
		     {
		       if (fu[i] > fb[i]) 
			 {
			   assert(u[i] >= lim_inf && u[i] <= lim_sup);
			   cx[i]=u[i];
			   fc[i]=fu[i];
			   converged[i] = TRUE;			  
			 }
		       else
			 {		   
			   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
			   param[i] = u[i];
			   if(param[i] > lim_sup) {param[i] = u[i] = lim_sup;}
			   if(param[i] < lim_inf) {param[i] = u[i] = lim_inf;}	  
			   state[i] = 1;		 
			 }		  
		     }
		   break;
		 case 2: 
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {		     
		       SHFT(bx[i],cx[i],u[i], cx[i]+MNBRAK_GOLD*(cx[i]-bx[i]));
		       state[i] = 2;
		     }	   
		   else
		     {
		       state[i] = 0;
		       SHFT(ax[i],bx[i],cx[i],u[i]);
		       SHFT(fa[i],fb[i],fc[i],fu[i]);
		     }
		   break;	   
		 case 3:		  
		   SHFT(fb[i],fc[i],fu[i], temp[i]);
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 default:
		   assert(0);
		 }
	     }
	 }
    }
   

   assert(0);
   rax_free(converged);
   rax_free(ulim);
   rax_free(u);
   rax_free(r);
   rax_free(q);
   rax_free(fu);
   rax_free(dum); 
   rax_free(temp);
   rax_free(state);   
   rax_free(endState);

  

   return(0);
}


/**********************************************************************************************************/
/* ALPHA PARAM ********************************************************************************************/




/* function for optimizing alpha parameter */


static void optAlpha(tree *tr, partitionList *pr, double modelEpsilon, linkageList *ll)
{
  int 
    i, 
    k,
    numberOfModels = ll->entries;
  
  double 
    lim_inf     = ALPHA_MIN,
    lim_sup     = ALPHA_MAX;
  double
    *startLH    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *startAlpha = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *endAlpha   = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_a         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_b         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_c         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fa        = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fb        = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fc        = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_param     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *result     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_x         = (double *)rax_malloc(sizeof(double) * numberOfModels);   

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
   int revertModel = 0;
#endif   

   evaluateGeneric(tr, pr, tr->start, TRUE);
   
   /* 
     at this point here every worker has the traversal data it needs for the 
     search, so we won't re-distribute it he he :-)
  */

  for(i = 0; i < numberOfModels; i++)
    {
      /* make sure that we want to optimize alpha here */

      assert(ll->ld[i].valid);
      
      /* get the representative alpha value for all partitions that share it.
	 In RAxML we have hard-coded all alpha parameters to be estimated separately.
	 hence linking alpha parameters has not been tested for a long long time, be cautious 
      */
      
      startAlpha[i] = pr->partitionData[ll->ld[i].partitionList[0]]->alpha;
      _a[i] = startAlpha[i] + 0.1;
      _b[i] = startAlpha[i] - 0.1;      
      if(_b[i] < lim_inf) 
	_b[i] = lim_inf;

      startLH[i] = 0.0;
      
      for(k = 0; k < ll->ld[i].partitions; k++)	
	{
	  /* sum over all per-partition log likelihood scores that share the same alpha parameter */
	  startLH[i] += pr->partitionData[ll->ld[i].partitionList[k]]->partitionLH;
	  
	  /* make sure that all copies of the linked alpha parameter have been updated appropriately */
	  assert(pr->partitionData[ll->ld[i].partitionList[0]]->alpha ==  pr->partitionData[ll->ld[i].partitionList[k]]->alpha);
	}
    }					  

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, -1, ALPHA_F, tr, pr, ll);
  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, ALPHA_F, -1, tr, pr, ll, lim_inf, lim_sup);

  for(i = 0; i < numberOfModels; i++)
    endAlpha[i] = result[i];
  
  /* if for some reason we couldn't improve the likelihood further, restore the old alpha value for all partitions sharing it */

  for(i = 0; i < numberOfModels; i++)
    {
      if(startLH[i] > endAlpha[i])
	{    	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {	      
	      pr->partitionData[ll->ld[i].partitionList[k]]->alpha = startAlpha[i];

	      makeGammaCats(pr->partitionData[ll->ld[i].partitionList[k]]->alpha, pr->partitionData[ll->ld[i].partitionList[k]]->gammaRates, 4, tr->useMedian);
	    }
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  revertModel++;
#endif
	}  
    }

  /* broadcast new alpha value to all parallel threads/processes */

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  if(revertModel > 0)
    masterBarrier(THREAD_COPY_ALPHA, tr, pr);
#endif
  
  rax_free(startLH);
  rax_free(startAlpha);
  rax_free(endAlpha);
  rax_free(result);
  rax_free(_a);
  rax_free(_b);
  rax_free(_c);
  rax_free(_fa);
  rax_free(_fb);
  rax_free(_fc);
  rax_free(_param);
  rax_free(_x);  

}

/* optimize rates in the Q matrix */

static void optRates(tree *tr, partitionList *pr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int 
    i, 
    k, 
    j, 
    pos,
    numberOfRates = ((states * states - states) / 2) - 1;
    
  double lim_inf = RATE_MIN;
  double lim_sup = RATE_MAX;
  double 
    *startRates;
  double 
    *startLH= (double *)rax_malloc(sizeof(double) * numberOfModels),
    *endLH  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_a     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_b     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_c     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fa    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fb    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fc    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_param = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *result = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_x     = (double *)rax_malloc(sizeof(double) * numberOfModels); 

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
   int revertModel = 0;
#endif

   assert(states != -1);

  startRates = (double *)rax_malloc(sizeof(double) * numberOfRates * numberOfModels);

  evaluateGeneric(tr, pr, tr->start, TRUE);
  
  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  /* get the initial rates */

  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  endLH[pos] = unlikely;
	  startLH[pos] = 0.0;

	  /* loop over all partitions that are linked via the corresponding GTR matrx, i.e., that share the same GTR matrix */

	  for(j = 0; j < ll->ld[i].partitions; j++)
	    {
	      int 
		index = ll->ld[i].partitionList[j];
	      
	      startLH[pos] += pr->partitionData[index]->partitionLH;
	      
	      for(k = 0; k < numberOfRates; k++)
		startRates[pos * numberOfRates + k] = pr->partitionData[index]->substRates[k];
	    }
	  pos++;
	}
    }  

  assert(pos == numberOfModels);
  
  /* 
     Now loop over all rates in the matrix, e.g., 5 if it's DNA, and 189 if it's protein data.
     Note that we do this on a rate by rate basis, i.e., first try to optimize a->c, a->g, etc.
   */

  for(i = 0; i < numberOfRates; i++)
    {     
      for(k = 0, pos = 0; k < ll->entries; k++)
	{
	  if(ll->ld[k].valid)
	    {
	      int index = ll->ld[k].partitionList[0];
	      _a[pos] = pr->partitionData[index]->substRates[i] + 0.1;
	      _b[pos] = pr->partitionData[index]->substRates[i] - 0.1;
	      
	      if(_a[pos] < lim_inf) _a[pos] = lim_inf;
	      if(_a[pos] > lim_sup) _a[pos] = lim_sup;
	      
	      if(_b[pos] < lim_inf) _b[pos] = lim_inf;
	      if(_b[pos] > lim_sup) _b[pos] = lim_sup;    
	      pos++;
	    }
	}                       	     

      assert(pos == numberOfModels);

      brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, i, RATE_F, tr, pr, ll);
      
      for(k = 0; k < numberOfModels; k++)
	{
	  assert(_a[k] >= lim_inf && _a[k] <= lim_sup);
	  assert(_b[k] >= lim_inf && _b[k] <= lim_sup);	  
	  assert(_c[k] >= lim_inf && _c[k] <= lim_sup);	    
	}      

      brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, RATE_F, i, tr, pr, ll, lim_inf, lim_sup);
	
      for(k = 0; k < numberOfModels; k++)
	endLH[k] = result[k];	           
      
      /* if the proposed new rate does not improve the likelihood undo the change */

      for(k = 0, pos = 0; k < ll->entries; k++)
	{
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  revertModel = 0;
#endif
	  if(ll->ld[k].valid)
	    { 
	      if(startLH[pos] > endLH[pos])
		{		  
		  for(j = 0; j < ll->ld[k].partitions; j++)
		    {
		      int index = ll->ld[k].partitionList[j];
		      pr->partitionData[index]->substRates[i] = startRates[pos * numberOfRates + i];

		      initReversibleGTR(tr, pr, index);
		    }
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
		  revertModel++;
#endif
		}
	      pos++;
	    }
	}

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      if(revertModel > 0)
	masterBarrier(THREAD_COPY_RATES, tr, pr);
#endif
      
      assert(pos == numberOfModels);
    }

 
  rax_free(startLH);
  rax_free(endLH);
  rax_free(result);
  rax_free(_a);
  rax_free(_b);
  rax_free(_c);
  rax_free(_fa);
  rax_free(_fb);
  rax_free(_fc);
  rax_free(_param);
  rax_free(_x);  
  rax_free(startRates);
}


/* figure out if all AA models have been assigned a joint GTR matrix */

static boolean AAisGTR(partitionList *pr)
{
  int i, count = 0;

  for(i = 0; i < pr->numberOfPartitions; i++)
    {
      if(pr->partitionData[i]->dataType == AA_DATA)
	{
	  count++;
	  if(pr->partitionData[i]->protModels != GTR)
	    return FALSE;
	}
    }

  if(count == 0)
    return FALSE;

  return TRUE;
}


/* generic substitiution matrix (Q matrix) optimization */

static void optRatesGeneric(tree *tr, partitionList *pr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    dnaPartitions = 0,
    aaPartitions  = 0,
    secondaryPartitions = 0,
    secondaryModel = -1,
    states = -1;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* 
     first optimize all rates in DNA data partition matrices. That's where we use the valid field in the 
     linkage list data structure. 
   */

  for(i = 0; i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
	{
	case DNA_DATA:	
	  states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	  ll->ld[i].valid = TRUE;
	  dnaPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  /* if we have dna partitions in our dataset, let's optimize all 5 rates in their substitution matrices */

  if(dnaPartitions > 0)
    optRates(tr, pr, modelEpsilon, ll, dnaPartitions, states);
  

  /* now we do the same thing for secondary structure models  */

   for(i = 0; i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
	{
	  /* 
	     we only have one type of secondary data models in one analysis, RAxML allows for only one secondary data 
	     partition !
	   */
	case SECONDARY_DATA_6:
	  states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	  secondaryModel = SECONDARY_DATA_6;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA_7: 
	  states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	  secondaryModel = SECONDARY_DATA_7;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA:
	  states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	  secondaryModel = SECONDARY_DATA;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:	
	case DNA_DATA:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }

   /* if there are secondary struct data partitions (and there can only be at most one), let's optimize 
      their rates, tehre are models with 6, 7 and 16 states */
   
   if(secondaryPartitions > 0)
     {
       assert(secondaryPartitions == 1);

       switch(secondaryModel)
	 {
	 case SECONDARY_DATA:
	   optRates(tr, pr, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_6:
	   optRates(tr, pr, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_7:
	   optRates(tr, pr, modelEpsilon, ll, secondaryPartitions, states);
	   break; 
	 default:
	   assert(0);
	 }
     }

  /* then AA for GTR */

   /* now if all AA partitions share a joint GTR subst matrix, let's do a joint estimate 
      of the 189 rates across all of them. Otherwise we don't need to optimize anything since 
      we will be using one of the fixed models like WAG, JTT, etc */

  if(AAisGTR(pr))
    {
      for(i = 0; i < ll->entries; i++)
	{
	  switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
	    {
	    case AA_DATA:
	      states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	      ll->ld[i].valid = TRUE;
	      aaPartitions++;
	      break;
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	      ll->ld[i].valid = FALSE;
	      break;
	    default:
	      assert(0);
	    }	 
	}

      assert(aaPartitions == 1);     
      
      optRates(tr, pr, modelEpsilon, ll, aaPartitions, states);
    }
  
  /* then multi-state data partitions

     now here we have to be careful, because every multi-state partition can actually 
     have a distinct number of states, so we will go to every multi-state partition separately,
     parallel efficiency for this will suck, but what the hell .....
  */

  if(tr->multiStateModel == GTR_MULTI_STATE)
    {     
      for(i = 0; i < ll->entries; i++)
	{
	  switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
	    {
	    case GENERIC_32:
	      {
		int k;
		
		states = pr->partitionData[ll->ld[i].partitionList[0]]->states;

		ll->ld[i].valid = TRUE;
		
		for(k = 0; k < ll->entries; k++)
		  if(k != i)
		    ll->ld[k].valid = FALSE;
		
		optRates(tr, pr, modelEpsilon, ll, 1, states);
	      }
	      break;
	    case AA_DATA:	    
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	    case GENERIC_64:
	      break;
	    default:
	      assert(0);
	    }	 
	}           
    }

  /* done with all partitions, so we can set all entries in the linkage list to valid again :-) */

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}





/*********************FUNCTIONS FOR PSR/CAT model of rate heterogeneity ***************************************/






static int catCompare(const void *p1, const void *p2)
{
 rateCategorize *rc1 = (rateCategorize *)p1;
 rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


static void categorizePartition(tree *tr, partitionList *pr, rateCategorize *rc, int model, int lower, int upper)
{
  int
    zeroCounter,
    i, 
    k;
  
  double 
    diff, 
    min;

  for (i = lower, zeroCounter = 0; i < upper; i++, zeroCounter++) 
      {
	double
	  temp = tr->patrat[i];

	int
	  found = 0;
	
	for(k = 0; k < pr->partitionData[model]->numberOfCategories; k++)
	  {
	    if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	      {
		found = 1;
		tr->rateCategory[i] = k; 
		break;
	      }
	  }
	
	if(!found)
	  {
	    min = fabs(temp - rc[0].rate);
	    tr->rateCategory[i] = 0;

	    for(k = 1; k < pr->partitionData[model]->numberOfCategories; k++)
	      {
		diff = fabs(temp - rc[k].rate);

		if(diff < min)
		  {
		    min = diff;
		    tr->rateCategory[i] = k;
		  }
	      }
	  }
      }

  for(k = 0; k < pr->partitionData[model]->numberOfCategories; k++)
    pr->partitionData[model]->perSiteRates[k] = rc[k].rate;
}


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

void optRateCatPthreads(tree *tr, partitionList *pr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid)
{
  int 
    model, 
    i;

  for(model = 0; model < pr->numberOfPartitions; model++)
    {      
      int 
	localIndex = 0;

      boolean 
	execute = ((tr->manyPartitions && isThisMyPartition(pr, tid, model)) || (!tr->manyPartitions));

      if(execute)
	for(i = pr->partitionData[model]->lower;  i < pr->partitionData[model]->upper; i++)
	  {
	    if(tr->manyPartitions || (i % n == tid))
	      {
	      
		double initialRate, initialLikelihood, 
		  leftLH, rightLH, leftRate, rightRate, v;
		const double epsilon = 0.00001;
		int k;	      
		
		tr->patrat[i] = tr->patratStored[i];     
		initialRate = tr->patrat[i];
		
		initialLikelihood = evaluatePartialGeneric(tr, pr, localIndex, initialRate, model); /* i is real i ??? */
		
		
		leftLH = rightLH = initialLikelihood;
		leftRate = rightRate = initialRate;
		
		k = 1;
		
		while((initialRate - k * lower_spacing > 0.0001) && 
		      ((v = evaluatePartialGeneric(tr, pr, localIndex, initialRate - k * lower_spacing, model))
		       > leftLH) && 
		      (fabs(leftLH - v) > epsilon))  
		  {	  
#ifndef WIN32
		    if(isnan(v))
		      assert(0);
#endif
		    
		    leftLH = v;
		    leftRate = initialRate - k * lower_spacing;
		    k++;	  
		  }      
		
		k = 1;
		
		while(((v = evaluatePartialGeneric(tr, pr, localIndex, initialRate + k * upper_spacing, model)) > rightLH) &&
		      (fabs(rightLH - v) > epsilon))    	
		  {
#ifndef WIN32
		    if(isnan(v))
		      assert(0);
#endif     
		    rightLH = v;
		    rightRate = initialRate + k * upper_spacing;	 
		    k++;
		  }           
		
		if(rightLH > initialLikelihood || leftLH > initialLikelihood)
		  {
		    if(rightLH > leftLH)	    
		      {	     
			tr->patrat[i] = rightRate;
			lhs[i] = rightLH;
		      }
		    else
		      {	      
			tr->patrat[i] = leftRate;
			lhs[i] = leftLH;
		      }
		  }
		else
		  lhs[i] = initialLikelihood;
		
		tr->patratStored[i] = tr->patrat[i];
		localIndex++;
	      }
	  }
      assert(localIndex == pr->partitionData[model]->width);
    }
}



#else


static void optRateCatModel(tree *tr, partitionList *pr, int model, double lower_spacing, double upper_spacing, double *lhs)
{
  int lower = pr->partitionData[model]->lower;
  int upper = pr->partitionData[model]->upper;
  int i;
  for(i = lower; i < upper; i++)
    {
      double initialRate, initialLikelihood, 
	leftLH, rightLH, leftRate, rightRate, v;
      const double epsilon = 0.00001;
      int k;
      
      tr->patrat[i] = tr->patratStored[i];     
      initialRate = tr->patrat[i];
      
      initialLikelihood = evaluatePartialGeneric(tr, pr, i, initialRate, model);
      
      
      leftLH = rightLH = initialLikelihood;
      leftRate = rightRate = initialRate;
      
      k = 1;
      
      while((initialRate - k * lower_spacing > 0.0001) && 
	    ((v = evaluatePartialGeneric(tr, pr, i, initialRate - k * lower_spacing, model))
	     > leftLH) && 
	    (fabs(leftLH - v) > epsilon))  
	{	  
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif
	  
	  leftLH = v;
	  leftRate = initialRate - k * lower_spacing;
	  k++;	  
	}      
      
      k = 1;
      
      while(((v = evaluatePartialGeneric(tr, pr, i, initialRate + k * upper_spacing, model)) > rightLH) &&
	    (fabs(rightLH - v) > epsilon))    	
	{
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif     
	  rightLH = v;
	  rightRate = initialRate + k * upper_spacing;	 
	  k++;
	}           
  
      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
	{
	  if(rightLH > leftLH)	    
	    {	     
	      tr->patrat[i] = rightRate;
	      lhs[i] = rightLH;
	    }
	  else
	    {	      
	      tr->patrat[i] = leftRate;
	      lhs[i] = leftLH;
	    }
	}
      else
	lhs[i] = initialLikelihood;
      
      tr->patratStored[i] = tr->patrat[i];
    }

}


#endif



/* 
   set scaleRates to FALSE everywhere such that 
   per-site rates are not scaled to obtain an overall mean rate 
   of 1.0
*/

void updatePerSiteRates(tree *tr, partitionList *pr, boolean scaleRates)
{
  int 
    i,
    model;

  if(pr->perGeneBranchLengths && pr->numberOfPartitions > 1)
    {            
      for(model = 0; model < pr->numberOfPartitions; model++)
	{
	  int 	       
	    width = pr->partitionData[model]->width,
	    lower = pr->partitionData[model]->lower,
	    upper = pr->partitionData[model]->upper;
	  
	  
	  if(scaleRates)
	    {
	      double 
		scaler = 0.0,       
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = 0; i < width; i++)
		{
		  int 
		    /* w = tr->aliaswgt[i]; */
		    w = pr->partitionData[model]->wgt[i]; 
		  
		  double
		    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	  
	      scaler = 1.0 / ((double)accRat);
	  	  
	      for(i = 0; i < pr->partitionData[model]->numberOfCategories; i++)
		pr->partitionData[model]->perSiteRates[i] *= scaler;

	      accRat = 0.0;	 
	      
	      for(i = 0; i < width; i++)
		{
		  int 
		    /* w = tr->[i]; */
		    w = pr->partitionData[model]->wgt[i]; 
		  
		  double
		    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}	         

	      accRat /= ((double)accWgt);	  

	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }
	  else
	    {
	      double 		   
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = 0; i < width; i++)
		{
		  int 
		    /* w = tr->aliaswgt[i]; */
		    w = pr->partitionData[model]->wgt[i]; 
		  
		  double
		    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	      
	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }

	  
#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  {
	    int 
	      localCount = 0;
	    
	    for(i = lower, localCount = 0; i < upper; i++, localCount++)
	      {	    	      
		pr->partitionData[model]->rateCategory[localCount] = tr->rateCategory[i];
	      }
	  }
#endif
	}
    }
  else
    {
      /* assert(0);  */
      int
	accWgt = 0;

      double 
	scaler = 0.0,       
	accRat = 0.0; 

      if(scaleRates)
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < pr->numberOfPartitions; model++)
	    {
	      int 
		localCount = 0,
		width = pr->partitionData[model]->width,
		lower = pr->partitionData[model]->lower,
		upper = pr->partitionData[model]->upper;
	      
	      for(i = lower ; i < upper; i++)
		{
		  int 
		    /* w = tr->aliaswgt[i]; */
		    w = pr->partitionData[model]->wgt[i]; 
		  
		  double
		    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		  ++localCount; 
		}
	    }
	  
	  accRat /= ((double)accWgt);
	  
	  scaler = 1.0 / ((double)accRat);
	  
	  for(model = 0; model < pr->numberOfPartitions; model++)
	    {
	      for(i = 0; i < pr->partitionData[model]->numberOfCategories; i++)
		pr->partitionData[model]->perSiteRates[i] *= scaler;
	    }

	  for(model = 0, accRat = 0.0; model < pr->numberOfPartitions; model++)
	    {
	      int 
		localCount = 0,
		lower = pr->partitionData[model]->lower,
		upper = pr->partitionData[model]->upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}
	    }           

	  accRat /= ((double)accWgt);	  

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
      else
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < pr->numberOfPartitions; model++)
	    {
	      int 
		localCount = 0,
		width = pr->partitionData[model]->width,
		lower = pr->partitionData[model]->lower,
		upper = pr->partitionData[model]->upper;
	      
	      for(i = 0; i < width; i++)
		{
		  int 
		    /* w = tr->aliaswgt[i]; */
		    w = pr->partitionData[model]->wgt[i]; 
		  
		  double
		    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		  ++localCount; 
		}
	    }
	  
	  accRat /=  (double)accWgt;

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
         
         /*
       for(model = 0; model < pr->numberOfPartitions; model++)
	{
	  int 
	    localCount = 0,
	    lower = pr->partitionData[model]->lower,
	    upper = pr->partitionData[model]->upper;

	}  */       
#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      for(model = 0; model < pr->numberOfPartitions; model++)
	{   	  	  	 
	  int 
	    localCount,
	    lower = pr->partitionData[model]->lower,
	    upper = pr->partitionData[model]->upper;
	  
	  for(i = lower, localCount = 0; i < upper; i++, localCount++)
	      pr->partitionData[model]->rateCategory[localCount] = tr->rateCategory[i];
	}
#endif
    }
  
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  masterBarrier(THREAD_COPY_RATE_CATS, tr, pr);
#endif               
}

static void optimizeRateCategories(tree *tr, partitionList *pr, int _maxCategories)
{
  assert(_maxCategories > 0);

  if(_maxCategories > 1)
    {
      double  
	temp,  
	lower_spacing, 
	upper_spacing,
	initialLH = tr->likelihood,	
	*ratStored = (double *)rax_malloc(sizeof(double) * tr->originalCrunchedLength),
	/**lhs =       (double *)malloc(sizeof(double) * tr->originalCrunchedLength),*/
	**oldCategorizedRates = (double **)rax_malloc(sizeof(double *) * pr->numberOfPartitions);

      int  
	i,
	k,
	maxCategories = _maxCategories,
	*oldCategory =  (int *)rax_malloc(sizeof(int) * tr->originalCrunchedLength),
	model,
	*oldNumbers = (int *)rax_malloc(sizeof(int) * pr->numberOfPartitions);
  
      assert(isTip(tr->start->number, tr->mxtips));         
      
      evaluateGeneric(tr, pr, tr->start, TRUE);

      if(tr->optimizeRateCategoryInvocations == 1)
	{
	  lower_spacing = 0.5 / ((double)(tr->optimizeRateCategoryInvocations));
	  upper_spacing = 1.0 / ((double)(tr->optimizeRateCategoryInvocations));
	}
      else
	{
	  lower_spacing = 0.05 / ((double)(tr->optimizeRateCategoryInvocations));
	  upper_spacing = 0.1 / ((double)(tr->optimizeRateCategoryInvocations));
	}
      
      if(lower_spacing < 0.001)
	lower_spacing = 0.001;
      
      if(upper_spacing < 0.001)
	upper_spacing = 0.001;
      
      tr->optimizeRateCategoryInvocations = tr->optimizeRateCategoryInvocations + 1;

      memcpy(oldCategory, tr->rateCategory, sizeof(int) * tr->originalCrunchedLength);	     
      memcpy(ratStored,   tr->patratStored, sizeof(double) * tr->originalCrunchedLength);

      for(model = 0; model < pr->numberOfPartitions; model++)
	{
	  oldNumbers[model]          = pr->partitionData[model]->numberOfCategories;

	  oldCategorizedRates[model] = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
	  
	  memcpy(oldCategorizedRates[model], pr->partitionData[model]->perSiteRates, tr->maxCategories * sizeof(double));
	}      
      
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      /*tr->lhs = lhs;*/
      tr->lower_spacing = lower_spacing;
      tr->upper_spacing = upper_spacing;
      masterBarrier(THREAD_RATE_CATS, tr, pr);
#else      
      for(model = 0; model < pr->numberOfPartitions; model++)
	optRateCatModel(tr, pr, model, lower_spacing, upper_spacing, tr->lhs);
#endif     

      for(model = 0; model < pr->numberOfPartitions; model++)
	{     
	  int 
	    where = 1,
	    found = 0,
	    width = pr->partitionData[model]->upper -  pr->partitionData[model]->lower,
	    upper = pr->partitionData[model]->upper,
	    lower = pr->partitionData[model]->lower;
	    
	  rateCategorize 
	    *rc = (rateCategorize *)rax_malloc(sizeof(rateCategorize) * width);		 
	
	  for (i = 0; i < width; i++)
	    {
	      rc[i].accumulatedSiteLikelihood = 0.0;
	      rc[i].rate = 0.0;
	    }  
	
	  rc[0].accumulatedSiteLikelihood = tr->lhs[lower];
	  rc[0].rate = tr->patrat[lower];
	
	  tr->rateCategory[lower] = 0;
	
	  for (i = lower + 1; i < upper; i++) 
	    {
	      temp = tr->patrat[i];
	      found = 0;
	    
	      for(k = 0; k < where; k++)
		{
		  if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
		    {
		      found = 1;						
		      rc[k].accumulatedSiteLikelihood += tr->lhs[i];	
		      break;
		    }
		}
	    
	      if(!found)
		{	    
		  rc[where].rate = temp;	    
		  rc[where].accumulatedSiteLikelihood += tr->lhs[i];	    
		  where++;
		}
	    }
	
	  qsort(rc, where, sizeof(rateCategorize), catCompare);
	
	  if(where < maxCategories)
	    {
	      pr->partitionData[model]->numberOfCategories = where;
	      categorizePartition(tr, pr, rc, model, lower, upper);
	    }
	  else
	    {
	      pr->partitionData[model]->numberOfCategories = maxCategories;
	      categorizePartition(tr, pr, rc, model, lower, upper);
	    }
	
	  rax_free(rc);
	}
      
      updatePerSiteRates(tr, pr, TRUE);

      evaluateGeneric(tr, pr, tr->start, TRUE);
      
      if(tr->likelihood < initialLH)
	{	 		  
	  for(model = 0; model < pr->numberOfPartitions; model++)
	    {
	      pr->partitionData[model]->numberOfCategories = oldNumbers[model];
	      memcpy(pr->partitionData[model]->perSiteRates, oldCategorizedRates[model], tr->maxCategories * sizeof(double));
	    }	      
	  
	  memcpy(tr->patratStored, ratStored, sizeof(double) * tr->originalCrunchedLength);
	  memcpy(tr->rateCategory, oldCategory, sizeof(int) * tr->originalCrunchedLength);	     
	  
	  updatePerSiteRates(tr, pr, FALSE);
	  
	  evaluateGeneric(tr, pr, tr->start, TRUE);

	  /* printf("REVERT: %1.40f %1.40f\n", initialLH, tr->likelihood); */

	  assert(initialLH == tr->likelihood);
	}
          
      for(model = 0; model < pr->numberOfPartitions; model++)
	rax_free(oldCategorizedRates[model]);
                   
      rax_free(oldCategorizedRates);
      rax_free(oldCategory);
      rax_free(ratStored);       
      /*     rax_free(lhs); */
      rax_free(oldNumbers);
    }
}
  

/************************* end of functions for CAT model of rate heterogeneity */




/*****************************************************************************************************/

/* reset all branche lengths in tree to default values */

void resetBranches(tree *tr)
{
  nodeptr  p, q;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
      for(i = 0; i < NUM_BRANCHES; i++)
	p->z[i] = defaultz;
	
      q = p->next;
      while(q != p)
	{	
	  for(i = 0; i < NUM_BRANCHES; i++)
	    q->z[i] = defaultz;	    
	  q = q->next;
	}
      p++;
    }
}

/* print the protein GTR substitution matrix to file. This is only printed if we are actually estimating a 
   GTR model for all protein data partitions in the data */

static void printAAmatrix(partitionList *pr, double epsilon)
{
  if(AAisGTR(pr))
    {
      int model;
      
      for(model = 0; model < pr->numberOfPartitions; model++)
	{
	  if(pr->partitionData[model]->dataType == AA_DATA)
	    {
	      char gtrFileName[1024];
	      char epsilonStr[1024];
	      FILE *gtrFile;
	      double *rates = pr->partitionData[model]->substRates;
	      double *f     = pr->partitionData[model]->frequencies;
	      double q[20][20];
	      int    r = 0;
	      int i, j;

	      assert(pr->partitionData[model]->protModels == GTR);

	      sprintf(epsilonStr, "%f", epsilon);

	      strcpy(gtrFileName, workdir);
	      strcat(gtrFileName, "RAxML_proteinGTRmodel.");
	      strcat(gtrFileName, run_id);
	      strcat(gtrFileName, "_");
	      strcat(gtrFileName, epsilonStr);

	      gtrFile = myfopen(gtrFileName, "wb");

	      for(i = 0; i < 20; i++)
		for(j = 0; j < 20; j++)
		  q[i][j] = 0.0;

	      for(i = 0; i < 19; i++)
		for(j = i + 1; j < 20; j++)
		  q[i][j] = rates[r++];

	      for(i = 0; i < 20; i++)
		for(j = 0; j <= i; j++)
		  {
		    if(i == j)
		      q[i][j] = 0.0;
		    else
		      q[i][j] = q[j][i];
		  }
	   
	      for(i = 0; i < 20; i++)
		{
		  for(j = 0; j < 20; j++)		
		    fprintf(gtrFile, "%1.80f ", q[i][j]);
		
		  fprintf(gtrFile, "\n");
		}
	      for(i = 0; i < 20; i++)
		fprintf(gtrFile, "%1.80f ", f[i]);
	      fprintf(gtrFile, "\n");

	      fclose(gtrFile);

	      printBothOpen("\nPrinted intermediate AA substitution matrix to file %s\n\n", gtrFileName);
	      
	      break;
	    }

	}	  
    }
}

/* 
   automatically compute the best protein substitution model for the dataset at hand.
 */

static void autoProtein(tree *tr, partitionList *pr)
{
  int 
    countAutos = 0,
    i,
    model;
    
   for(model = 0; model < pr->numberOfPartitions; model++)
     if(pr->partitionData[model]->protModels == AUTO)
       countAutos++;
  
  if(countAutos > 0)
    {
      int 
	numProteinModels = AUTO,
	*bestIndex = (int*)rax_malloc(sizeof(int) * pr->numberOfPartitions);

      double
	*bestScores = (double*)rax_malloc(sizeof(double) * pr->numberOfPartitions);

      /*printf("Entry: %f\n", tr->likelihood);*/
      
      for(model = 0; model < pr->numberOfPartitions; model++)
	{
	  bestIndex[model] = -1;
	  bestScores[model] = unlikely;
	}

      for(i = 0; i < numProteinModels; i++)
	{
	  for(model = 0; model < pr->numberOfPartitions; model++)
	    {	   
	      if(pr->partitionData[model]->protModels == AUTO)
		{
		  pr->partitionData[model]->autoProtModels = i;

		  initReversibleGTR(tr, pr, model);
		}
	    }
	  
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  masterBarrier(THREAD_COPY_RATES, tr, pr);
#endif
	  
	  resetBranches(tr);
	  evaluateGeneric(tr, pr, tr->start, TRUE);
	  treeEvaluate(tr, pr, 16);// 0.5 * 32 = 16.0

	  for(model = 0; model < pr->numberOfPartitions; model++)
	    {
	      if(pr->partitionData[model]->protModels == AUTO)
		{
		  if(pr->partitionData[model]->partitionLH > bestScores[model])
		    {
		      bestScores[model] = pr->partitionData[model]->partitionLH;
		      bestIndex[model] = i;		      
		    }
		}

	    }       
	}

      printBothOpen("\n\n");
      
      for(model = 0; model < pr->numberOfPartitions; model++)
	{	   
	  if(pr->partitionData[model]->protModels == AUTO)
	    {
	      pr->partitionData[model]->autoProtModels = bestIndex[model];

	      initReversibleGTR(tr, pr, model);

	      printBothOpen("Partition: %d best-scoring AA model: %s likelihood %f\n", model, protModels[pr->partitionData[model]->autoProtModels], bestScores[model]);
	    }	 
	}
      
      printBothOpen("\n\n");
            
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      masterBarrier(THREAD_COPY_RATES, tr, pr);
#endif

      resetBranches(tr);
      evaluateGeneric(tr, pr, tr->start, TRUE);
      treeEvaluate(tr, pr, 16); // 0.5 * 32 = 16
      
      /*printf("Exit: %f\n", tr->likelihood);*/
      
      rax_free(bestIndex);
      rax_free(bestScores);
    }
}

/* iterative procedure for optimizing all model parameters */

void modOpt(tree *tr, partitionList *pr, double likelihoodEpsilon)
{ 
  int i, catOpt = 0; 
  double 
    currentLikelihood,
    modelEpsilon = 0.0001;

  /* linkage lists for alpha, p-invar has actually been ommitted in this version of the code 
     and the GTR subst matrices */

  linkageList *alphaList;
  linkageList *invarList;
  linkageList *rateList; 

  int *unlinked = (int *)rax_malloc(sizeof(int) * pr->numberOfPartitions);

  modelEpsilon = 0.0001;

  /* set the integer vector for linking parameters to all parameters being unlinked, 
     i.e. a separate/independent estimate of parameters will be conducted for eahc partition */

  for(i = 0; i < pr->numberOfPartitions; i++)
    unlinked[i] = i;

  /* alpha parameters and p-invar parameters are unlinked.
     this is the point where I actually hard-coded this in RAxML */

  alphaList = initLinkageList(unlinked, pr);
  invarList = initLinkageList(unlinked, pr);

  /* call the dedicated function for linking the GTR matrix across all AA data partitions 
     If we have only DNA data all GTR matrix estimates will be unlinked.
     */

  rateList  = initLinkageListGTR(pr);

  tr->start = tr->nodep[1];

  do
  {           
    //printBothOpen("cur LH: %f\n", tr->likelihood);
    currentLikelihood = tr->likelihood;     

    optRatesGeneric(tr, pr, modelEpsilon, rateList);

    evaluateGeneric(tr, pr, tr->start, TRUE);

    autoProtein(tr, pr);

    treeEvaluate(tr, pr, 2); // 0.0625 * 32 = 2.0

    switch(tr->rateHetModel)
    {
      case GAMMA:      
        optAlpha(tr, pr, modelEpsilon, alphaList);
        evaluateGeneric(tr, pr, tr->start, TRUE);
        treeEvaluate(tr, pr, 3); // 0.1 * 32 = 3.2
        break;
      case CAT:
        if(catOpt < 3)
        {	      	     	     
          optimizeRateCategories(tr, pr, tr->categories);
          catOpt++;
        }
        break;	  
      default:
        assert(0);
    }                   



    printAAmatrix(pr, fabs(currentLikelihood - tr->likelihood));
  }
  while(fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);  

  rax_free(unlinked);
  freeLinkageList(alphaList);
  freeLinkageList(rateList);
  freeLinkageList(invarList);  
}

