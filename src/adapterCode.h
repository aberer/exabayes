#ifndef _ADAPTER_CODE_H
#define _ADAPTER_CODE_H

#include "config.h"



#if HAVE_PLL == 1 
#define GET_NUM_PARTITIONS(tr) gAInfo.partitions.numberOfPartitions
#define GET_PARTITION(tr, model) (*(gAInfo.partitions.partitionData[model]))
#define HAS_PERGENE_BL(tr) (gAInfo.partitions.perGeneBranchLengths == TRUE)


				    /* NOTICE : here we assume, that in case of per-gene branch lengths, we have ONE BL for each model */
#define GET_NUM_BRANCHES(tr) ( gAInfo.partitions.perGeneBranchLengths ? gAInfo.partitions.numberOfPartitions : 1  ) 
#define GET_FRACCHANGE(tr,model) (gAInfo.partitions.partitionData[model]->fracchange)

#define GET_PCONTR(tr,model) (gAInfo.partitions.partitionData[model]->partitionContribution)
#define GET_PLH(tr,model) (gAInfo.partitions.partitionData[model]->partitionLH)

#define GET_EXECMODEL(tr,model) gAInfo.partitions.partitionData[model]->executeModel

#else 



#define GET_NUM_PARTITIONS(tr) tr->NumberOfModels
#define GET_PARTITION(tr,model) (tr->partitionData[model])
#define HAS_PERGENE_BL(tr) (tr->numBranches > 1)
#define GET_NUM_BRANCHES(tr) tr->numBranches
#define GET_FRACCHANGE(tr,model) tr->fracchanges[model]
#define GET_PCONTR(tr, model)  tr->partitionContributions[model] 
#define GET_PLH(tr,model) tr->perPartitionLH[model]
#define GET_EXECMODEL(tr,model) tr->executeModel[model]

#endif


void exa_newViewGeneric(tree *tr, nodeptr p, boolean masked); 
void exa_initReversibleGTR( tree *tr,int model); 
void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 


#endif
