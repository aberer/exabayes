#include "axml.h"

#include "proposalStructs.h"
#include "globals.h"

#include "main-common.h"


static void getxnodeBips (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->xBips || (s = s->next)->xBips)
    {
      p->xBips = s->xBips;
      s->xBips = 0;
    }

  assert(p->xBips);
}

static void newviewBipartitions(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength)
{
  
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q = p->next->back, 
      r = p->next->next->back;
    
    unsigned int       
      *vector = bitVectors[p->number],
      *left  = bitVectors[q->number],
      *right = bitVectors[r->number];
    unsigned 
      int i;      

    
    while(!p->xBips)
      {	
	if(!p->xBips)
	  getxnodeBips(p);
      }

    p->hash = q->hash ^ r->hash;

    if(isTip(q->number, numsp) && isTip(r->number, numsp))
      {		
	for(i = 0; i < vectorLength; i++)
	  vector[i] = left[i] | right[i];	  	
      }
    else
      {	
	if(isTip(q->number, numsp) || isTip(r->number, numsp))
	  {
	    if(isTip(r->number, numsp))
	      {	
		nodeptr tmp = r;
		r = q;
		q = tmp;
	      }	   
	    	    
	    while(!r->xBips)
	      {
		if(!r->xBips)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    	 
	  }
	else
	  {	    
	    while((!r->xBips) || (!q->xBips))
	      {
		if(!q->xBips)
		  newviewBipartitions(bitVectors, q, numsp, vectorLength);
		if(!r->xBips)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   	    	    	    	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	 
	  }

      }     
  }     
}


void insertAndCount(tree *tr, unsigned int *bitVector, hashtable *h, hashNumberType position, int chainId)
{     
  int vectorLength = ( tr->mxtips / 32  ) + ((tr->mxtips % 32) == 0 ?0 : 1) ; ; 
  
  
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    { 
	      e->treeVector[chainId]++; 
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 

      e->bitVector = (unsigned int*) exa_malloc_aligned(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));
      
      
      e->treeVector = exa_calloc(gAInfo.numberOfRuns, sizeof(nat)); 
      e->treeVector[chainId]++; 
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 

      e->bitVector = (unsigned int*)exa_malloc_aligned(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      e->treeVector = exa_calloc(gAInfo.numberOfRuns, sizeof(nat)); 
      e->treeVector[chainId]++; 

      memcpy(e->bitVector, bitVector, sizeof(nat) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}


void extractBipartitions(tree *tr, unsigned int **bitVectors, nodeptr p, hashtable *h, int *cnt, int chainId)
{
  int vectorLength = ( tr->mxtips / 32  ) + ((tr->mxtips % 32) == 0 ?0 : 1) ; ; 

  if(isTip(p->number, tr->mxtips))
    return;
  else
    {
      nodeptr 
	q = p->next; 

      do 
	{
	  extractBipartitions(tr, bitVectors, q->back, h, cnt, chainId); 
	  q = q->next;
	}
      while(q != p);
           
      newviewBipartitions(bitVectors, p, tr->mxtips, vectorLength);
      
      assert(p->xBips);

      if( NOT (isTip(p->back->number, tr->mxtips)))
	{
	  nat 
	    *toInsert  = bitVectors[p->number];

	  hashNumberType 
	    position = p->hash % h->tableSize;
	 
	  assert(!(toInsert[0] & 1));
	  *cnt += 1; 
	  insertAndCount(tr, toInsert, h, position, chainId); 
	} 
    }
}



void resetBitVectors(tree *tr)
{
  int bvLength = tr->mxtips / 32 +  ( (  tr->mxtips % 32  ) == 0 ? 0 : 1 )  ; 

  for(int i = 1 ; i < 2 * tr->mxtips ; ++i)
    {
      unsigned int *bv = tr->bitVectors[i]; 
      memset(bv, 0,  bvLength * sizeof(nat)); 
    }
}




void addBipartitionsToHash(tree *tr, state *chain)
{
  resetBitVectors(tr);

  int cnt = 0; 
  extractBipartitions(tr, tr->bitVectors, tr->start->back, chain->bvHash, &cnt, chain->id / gAInfo.numberCoupledChains);
  assert(cnt == tr->mxtips -3); 
}



/* TODO does not diagnose anything currently    */
boolean convergenceDiagnosticDummy(state *allChains, int numChains)
{
  boolean  result = TRUE; 
  
  for(int i = 0; i < numChains; ++i)
    {
      state *curChain = allChains + i; 
      if(curChain->currentGeneration < curChain->numGen)
	result = FALSE; 
    }
  
  return result; 
}





boolean averageDeviationOfSplitFrequencies(state *allChains)
{
  double asdsf = 0; 

  state *aChain = allChains + 0; 
  int numTaxa = aChain->tr->mxtips; 
  hashtable *ht = aChain->bvHash;

  int *numSampled = exa_calloc(gAInfo.numberOfRuns, sizeof(int)); 
  
  for(int i = 0; i < ht->tableSize; ++i)
    {
      entry *e = ht->table[i]; 
      
      while(e != NULL)
	{
	  for(int j = 0; j < gAInfo.numberOfRuns; ++j)
	    numSampled[j] += e->treeVector[j]; 

	  e = e->next; 
	}
    }
  

  /* asserting */
  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    assert(numSampled[0] == numSampled[i]); 

  assert(numSampled[0]  %  (numTaxa -3 )  == 0) ; 

  int treesSampled = numSampled[0] / (numTaxa - 3  ); 

  
  int cnt = 0; 
  int cntRelevant = 0; 
  for(int i = 0; i < ht->tableSize; ++i)
    {
      entry *e  = ht->table[i]; 

      while(e != NULL)
	{
	  /* check, if we can ignore this bipartition */
	  boolean relevant = FALSE; 
	  for(int j = 0; j < gAInfo.numberOfRuns; ++j)
	    if( e->treeVector[j] >  ((double)ASDSF_FREQ_CUTOFF * (double)treesSampled) )
	      relevant = TRUE; 

	  if(relevant)
	    {
	      double mu = 0, 
		sd = 0; 

	      for(int j = 0; j < gAInfo.numberOfRuns; ++j)
		mu += (double)e->treeVector[j] / (double)treesSampled; 

	      mu /= gAInfo.numberOfRuns; 

	      for(int j = 0; j < gAInfo.numberOfRuns; ++j)
		sd += pow((((double)e->treeVector[j]  / (double)treesSampled)  - mu ) ,2); 

#ifdef ASDSF_BE_VERBOSE
	      if(processID == 0)
		{
		  printf("RELEVANT: "); 
		  for(int  j = 0;  j < gAInfo.numberOfRuns; j++)
		    printf("%d,", e->treeVector[j]); 
		  printf("\tmu=%f\tsd=%f\n", mu, sd); 
		}
#endif

	      asdsf +=  sqrt(sd / gAInfo.numberOfRuns); 
	      cntRelevant++; 	    
	    }
	  else 
	    {
#ifdef ASDSF_BE_VERBOSE
	      if(processID == 0)
		{
		  printf("not relevant: "); 
		  for(int j = 0; j < gAInfo.numberOfRuns; ++j)
		    printf("%d,", e->treeVector[j]); 
		  printf("\n"); 
		}
#endif
	    }
	  

	  e = e->next; 
	  ++cnt; 
	}
    }

  assert(cnt == ht->entryCount); 
  asdsf /= cntRelevant; 

#ifdef ASDSF_BE_VERBOSE
  if(processID == 0)
    printf("CONVERGENCE: ASDSF = %f\n", asdsf); 
#endif

  exa_free(numSampled); 

  return (asdsf < ASDSF_CONVERGENCE_CRITERION ) ; 
}



boolean convergenceDiagnostic(state *allChains)
{
  if(gAInfo.numberOfRuns > 1)
    return averageDeviationOfSplitFrequencies(allChains); 
  else 
    return allChains[0].currentGeneration > allChains[0].numGen; 
}

