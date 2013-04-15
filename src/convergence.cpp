

#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "main-common.h"
#include "TreeAln.hpp" 


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

void printBv(unsigned int *bv, int bvlen)
{
  if(processID == 0)
    {
      for(int i = 0; i < bvlen; ++i)
	{
	  int pos = i / 32; 
	  int relPos = i % 32; 
	if( ( bv[pos ] & (1 << (relPos)) ) > 0 )
	  printf("1"); 
	else 
	  printf("0");
	}
    }
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
	  int i;
	  
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
      
      
      e->treeVector = (nat*)exa_calloc(gAInfo.numberOfRuns, sizeof(nat)); 
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

      e->treeVector = (nat*)exa_calloc(gAInfo.numberOfRuns, sizeof(nat)); 
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
  /* i think this function is not even necessary =/   */

  int bvLength = tr->mxtips / 32 +  ( (  tr->mxtips % 32  ) == 0 ? 0 : 1 )  ; 

  for(int i = 1 ; i < 2 * tr->mxtips ; ++i)
    {
      unsigned int *bv = tr->bitVectors[i]; 
      memset(bv, 0,  bvLength * sizeof(nat)); 
    }
  

  for(int i = 1 ; i <= tr->mxtips; ++i )
    {
      int pos = i / 32, 
	relPos = i % 32  ;
      tr->bitVectors[i][ pos  ] |=  (1 << relPos) ; 
    }
}




void addBipartitionsToHash(tree *tr, state *chain)
{
  resetBitVectors(tr);
  
  int cnt = 0; 
  extractBipartitions(tr, tr->bitVectors, tr->nodep[1]->back, gAInfo.bvHash, &cnt, chain->id / gAInfo.numberCoupledChains);
  assert(cnt == tr->mxtips -3); 
}



/* TODO does not diagnose anything currently    */
boolean convergenceDiagnosticDummy(state *allChains, int numChains)
{
  boolean  result = TRUE; 
  
  for(int i = 0; i < numChains; ++i)
    {
      state *curChain = allChains + i; 
      if(curChain->currentGeneration < gAInfo.numGen)
	result = FALSE; 
    }
  
  return result; 
}


boolean averageDeviationOfSplitFrequencies(state *allChains)
{
  double asdsf = 0; 

  state *aChain = allChains + 0; 
  int numTaxa = aChain->traln->getTr()->mxtips; 
  hashtable *ht = gAInfo.bvHash; 
  
  int *numSampled = (int*)exa_calloc(gAInfo.numberOfRuns, sizeof(int)); 
  
  entry** bvs = (entry**)exa_malloc(sizeof(entry*) * ht->entryCount); 
  nat ctr = 0; 
  
  for(nat i = 0; i < ht->tableSize; ++i)
    {
      entry *e = ht->table[i]; 
      
      while(e != NULL)
	{
	  for(int j = 0; j < gAInfo.numberOfRuns; ++j)
	    {
	      numSampled[j] += e->treeVector[j]; 
	    }
	  
	  bvs[ctr] = e;
	  ctr++; 
	  e = e->next; 
	}
    }
  assert(ctr == ht->entryCount); 

  /* asserting */
  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    assert(numSampled[0] == numSampled[i]); 
  assert(numSampled[0]  %  (numTaxa -3 )  == 0) ; 

  double treesSampled = numSampled[0] / (numTaxa - 3  ); 
#ifdef DEBUG_ASDSF_PRINT_ALL_BIPS
  if(processID == 0)
    printf("trees sampled: %g\n", treesSampled)  ;
#endif

  int cntRelevant = 0; 
  for(nat i = 0; i < ht->entryCount ; ++i)
    {
      entry *e = bvs[i]; 

      boolean relevant = FALSE; 
      for(int j = 0; j < gAInfo.numberOfRuns; ++j)
	if( e->treeVector[j] >  (gAInfo.asdsfIgnoreFreq * treesSampled) )
	  relevant = TRUE; 

	  if(relevant)
	    {
	      double 
		n = gAInfo.numberOfRuns,
		sd = 0; 

	      double mu = 0; 
	      for(int j = 0; j < n; ++j)
		mu += (double)e->treeVector[j] / treesSampled; 
	      mu /= n ; 

	      for(int j = 0; j < n; ++j)
		sd += pow((( ((double)e->treeVector[j]) / treesSampled) - mu ) ,2); 
	      sd = sqrt(sd/n) ; 

#ifdef DEBUG_ASDSF_PRINT_ALL_BIPS
	      if(processID == 0)
		{
		  printBv(e->bitVector, numTaxa); 
		  printf("\t");
		  for(int  j = 0;  j < n ; j++)
		    printf("%d,", e->treeVector[j]); 
		  printf("\tmu=%f\tsd=%f\n", mu, sd); 
		}
#endif 
	      asdsf += sd ; 
	      cntRelevant++; 	    
	    }
	  else 
	    {
#ifdef DEBUG_ASDSF_PRINT_ALL_BIPS
	      /* if(processID == 0) */
	      /* 	{ */
	      /* 	  printf("NOT\t");  */
	      /* 	  printBv(e->bitVector, numTaxa);  */
	      /* 	  printf("\t");  */
	      /* 	  for(int j = 0; j < gAInfo.numberOfRuns; ++j) */
	      /* 	    printf("%d,", e->treeVector[j]); */
	      /* 	  printf("\n"); */
	      /* 	} */
#endif
	    }
      
    }

  asdsf /= cntRelevant; 

  if(processID == 0)
    PRINT("CONVERGENCE: ASDSF = %f\n\n", asdsf); 

  exa_free(numSampled); 

  return (asdsf <  gAInfo.asdsfConvergence ) ; 
}



boolean convergenceDiagnostic(state *allChains)
{
  if(gAInfo.numberOfRuns > 1)
    return averageDeviationOfSplitFrequencies(allChains); 
  else 
    return allChains[0].currentGeneration > gAInfo.numGen; 
}

