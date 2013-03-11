
#include "common.h"
#include "axml.h"

#include "globals.h"

#include "config.h"

/** 
    Adapter methods that allow to either use examl or the pll for the
    build go in here
*/




/* int getNumberOfPartitions(tree *tr) */
/* { */
/* #if HAVE_PLL == 1  */
/*   assert(0);  */
/* #else  */
/*   return tr->NumberOfModels;  */
/* #endif */
/* } */


/* pInfo* getPartition(tree *tr, int model) */
/* { */
/* #if HAVE_PLL == 1  */
/*   assert(0);  */
/* #else  */
/*   return tr->partitionData + model;  */
/* #endif */
/* } */


void exa_newViewGeneric(tree *tr, nodeptr p, boolean masked)
{
#if HAVE_PLL == 1 
  newviewGeneric(tr, &(gAInfo.partitions),     p, masked); 
#else 
  newviewGeneric(tr, p, masked); 
#endif 
} 


void exa_initReversibleGTR( tree *tr,int model)
{
#if HAVE_PLL == 1
  initReversibleGTR(tr, &(gAInfo.partitions), model); 
#else 
  initReversibleGTR(tr, model); 
#endif
}


void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL == 1 
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}
