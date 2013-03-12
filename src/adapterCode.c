
#include "common.h"
#include "axml.h"

#include "globals.h"

#include "config.h"

/** 
    Adapter methods that allow to either use examl or the pll for the
    build go in here
*/







void setExecModel(tree *tr, int num,boolean value)
{
#if HAVE_PLL == 1 
  gAInfo.partitions.partitionData[num]->executeModel = value; 
#else 
  tr->executeModel[num] = value; 
#endif
}

boolean getExecModel(tree *tr, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.partitionData[num]->executeModel;   
#else 
  return tr->executeModel[num]; 
#endif
} 


void setPLH(tree *tr, int num, double value)
{
#if HAVE_PLL == 1 
  gAInfo.partitions.partitionData[num]->partitionLH = value; 
#else 
tr->perPartitionLH[num] = value; 
#endif
}


double getPLH(tree *tr, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.partitionData[num]->partitionLH; 
#else 
  return tr->perPartitionLH[num]; 
#endif
}



double getPcontr(tree *tr, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.partitionData[num]->partitionContribution; 
#else  
  return tr->partitionContributions[num]; 
#endif
}


double getFracChange(tree *tr, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.partitionData[num]->fracchange; 
#else
  return tr->fracchanges[num]; 
#endif
}


void setFracChange(tree *tr, int num, double value)
{
#if HAVE_PLL == 1 
  gAInfo.partitions.partitionData[num]->fracchange = value; 
#else 
  tr->fracchanges[num] = value; 
#endif
}



int getNumBranches(tree *tr)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.perGeneBranchLengths ? gAInfo.partitions.numberOfPartitions : 1; 
#else 
  return tr->numBranches; 
#endif
}



boolean hasPergeneBL(tree *tr)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.perGeneBranchLengths == TRUE; 
#else 
  return tr->numBranches > 1; 
#endif
}




int getNumberOfPartitions(tree *tr) 
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.numberOfPartitions; 
#else 
  return tr->NumberOfModels; 
#endif
}



pInfo* getPartition(tree *tr, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions.partitionData[num]; 
#else 
  return &(tr->partitionData[num]); 
#endif
} 





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



void exa_evaluateGeneric(tree *tr, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL == 1 
  evaluateGeneric(tr, &(gAInfo.partitions), start, fullTraversal); 
#else 
  evaluateGeneric(tr, start, fullTraversal); 
#endif  
}
