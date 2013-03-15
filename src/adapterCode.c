
#include "common.h"
#include "axml.h"

#include "proposalStructs.h"
#include "globals.h"

#include "config.h"

/** 
    Adapter methods that allow to either use examl or the pll for the
    build go in here
*/





void setNumberOfPartitions(tree *tr, int num)
{
#if HAVE_PLL == 1
  /* NOTICE if really something should be changed, do it! */   
#else 
  tr->NumberOfModels = num; 
#endif
}



void setNumbranches(tree *tr, int num)
{
#if HAVE_PLL == 1 
  /* NOTICE if really something should be changed, do it! */   
#else 
  tr->numBranches = num; 
#endif
}




void setExecModel(state *chain, int num,boolean value)
{
#if HAVE_PLL == 1 
  gAInfo.partitions[chain->couplingId].partitionData[num]->executeModel = value; 
#else 
  chain->tr->executeModel[num] = value; 
#endif
}

boolean getExecModel(state *chain, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[chain->couplingId].partitionData[num]->executeModel;   
#else 
  return chain->tr->executeModel[num]; 
#endif
} 


void setPLH(state *chain, int num, double value)
{
#if HAVE_PLL == 1 
  gAInfo.partitions[chain->couplingId].partitionData[num]->partitionLH = value; 
#else 
  chain->tr->perPartitionLH[num] = value; 
#endif
}


double getPLH(state *chain, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[chain->couplingId].partitionData[num]->partitionLH; 
#else 
  return chain->tr->perPartitionLH[num]; 
#endif
}



double getPcontr(state *chain, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[chain->couplingId].partitionData[num]->partitionContribution; 
#else  
  return chain->tr->partitionContributions[num]; 
#endif
}


double getFracChange(state *chain, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[chain->couplingId].partitionData[num]->fracchange; 
#else
  return chain->tr->fracchanges[num]; 
#endif
}


void setFracChange(state *chain, int num, double value)
{
#if HAVE_PLL == 1 
  gAInfo.partitions[chain->couplingId].partitionData[num]->fracchange = value; 
#else 
  chain->tr->fracchanges[num] = value; 
#endif
}



int getNumBranches(tree *tr)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[0].perGeneBranchLengths ? gAInfo.partitions[0].numberOfPartitions : 1; 
#else 
  return tr->numBranches; 
#endif
}



boolean hasPergeneBL(tree *tr)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[0].perGeneBranchLengths == TRUE; 
#else 
  return tr->numBranches > 1; 
#endif
}




int getNumberOfPartitions(tree *tr) 
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[0].numberOfPartitions; 
#else 
  return tr->NumberOfModels; 
#endif
}




pInfo* getPartition(state *chain, int num)
{
#if HAVE_PLL == 1 
  return gAInfo.partitions[chain->couplingId].partitionData[num]; 
#else 
  return chain->tr->partitionData + num; 
#endif
} 





void exa_newViewGeneric(state *chain, nodeptr p, boolean masked)
{
#if HAVE_PLL == 1 
  newviewGeneric(chain->tr, gAInfo.partitions + chain->couplingId,     p, masked); 
#else 
  newviewGeneric(chain->tr, p, masked); 
#endif 
} 


void exa_initReversibleGTR( state *chain,int model)
{
#if HAVE_PLL == 1
  initReversibleGTR(chain->tr, gAInfo.partitions + chain->couplingId, model); 
#else 
  initReversibleGTR(chain->tr, model); 
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



void exa_evaluateGeneric(state *chain, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL == 1 
  evaluateGeneric(chain->tr, gAInfo.partitions + chain->couplingId, start, fullTraversal); 
#else 
  evaluateGeneric(chain->tr, start, fullTraversal); 
#endif  
}


