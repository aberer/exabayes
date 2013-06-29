/**
   @file BipartitionHash.cpp

   @brief This file wraps good ol' RAxML bipartitionList.c into a
   bit-vector-based hash of bipartitions

   An important limitation of this hash is, that the mapping of taxa
   to numbers must always be the same 
 */

#include <cstring>
#include <cassert>

#include "BipartitionHash.hpp"
#include "TreeAln.hpp"
// #include "output.h"


/** 
    @brief allows to count bipartition frequencies

    @param numTax -- the number of taxa  in all trees that will hash their bipartitions into this object
    @param numRun -- the number of counters (usually different runs) for each bipartition 
 */ 
BipartitionHash::BipartitionHash(int _numTax, int numRuns)
  : numSlots(numRuns)
  , numTax(_numTax)
  , randomHash(numTax+1)
  , taxonNames(0)
{
  bitvectors = initBitVector(numTax, &vectorLength);
  h = initHashTable(10 * numTax);    
  initializeRandomHash(); 
}



BipartitionHash::~BipartitionHash()
{
  for(nat i = 1; i < 2 * numTax ;++i)
    exa_free(bitvectors[i]); 
  exa_free(bitvectors); 

  freeHashTable(h);
}



void BipartitionHash::printBv(nat *bv)
{
  int bvlen = vectorLength; 
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



void BipartitionHash::initializeRandomHash()
{
  // perfectly fine, we just want some deterministic randomness and it
  // is better not to use our actual RNG for that

  srand(1); 	
  
  for(nat i = 0; i <  numTax + 1 ; ++i)
    randomHash[i] = rand();
}




void BipartitionHash::extractBipartitions(TreeAln &traln, nodeptr p, int *cnt, int chainId)
{
  tree *tr = traln.getTr();

  if(isTip(p->number, tr->mxtips))
    return;
  else
    {
      nodeptr 
	q = p->next; 

      do 
	{
	  extractBipartitions(traln, q->back, cnt, chainId); 
	  q = q->next;
	}
      while(q != p);
           
      newviewBipartitions( p);

      if( NOT (isTip(p->back->number, tr->mxtips)))
	{
	  nat 
	    *toInsert  = bitvectors[p->number];

	  hashNumberType 
	    position = p->hash % h->tableSize;

	  assert(!(toInsert[0] & 1));

	  *cnt += 1; 
	  insertAndCount(traln, toInsert, position, chainId); 
	} 
    }
}



void BipartitionHash::resetBitVectors()
{
  for(nat i = 1 ; i < 2 * numTax ; ++i)
    {
      nat *bv = bitvectors[i]; 
      memset(bv, 0,  vectorLength * sizeof(nat)); 
    }

  for(nat i = 1 ; i <= numTax; ++i )
    {
      int pos = i / 32, 
	relPos = i % 32  ;      
      bitvectors[i][ pos  ] |=  (1 << relPos) ;       
    }
}




void BipartitionHash::getxnodeBips (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->xBips || (s = s->next)->xBips)
    {
      p->xBips = s->xBips;
      s->xBips = 0;
    }

  assert(p->xBips);
}



void BipartitionHash::insertAndCount(TreeAln &traln, nat *bitVector, hashNumberType position, int chainId)
{     

  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{
	  nat i ; 
	  for( i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if( i == vectorLength)
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
      
      
      e->treeVector = (nat*)exa_calloc(numSlots, sizeof(nat)); 
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

      e->treeVector = (nat*)exa_calloc(numSlots, sizeof(nat)); 
      e->treeVector[chainId]++; 

      memcpy(e->bitVector, bitVector, sizeof(nat) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}




void BipartitionHash::newviewBipartitions(nodeptr p)
{  
  if(isTip(p->number, numTax))
    return;

  nodeptr 
    q = p->next->back, 
    r = p->next->next->back;
    
  unsigned int       
    *vector = bitvectors[p->number],
    *left  = bitvectors[q->number],
    *right = bitvectors[r->number];
  unsigned 
    int i;      

    
  while(!p->xBips)
    {	
      if(!p->xBips)
	getxnodeBips(p);
    }

  p->hash = q->hash ^ r->hash;

  if(isTip(q->number, numTax) && isTip(r->number, numTax))
    {		
      for(i = 0; i < vectorLength; i++)
	vector[i] = left[i] | right[i];	  	
    }
  else
    {	
      if(isTip(q->number, numTax) || isTip(r->number, numTax))
	{
	  if(isTip(r->number, numTax))
	    {	
	      nodeptr tmp = r;
	      r = q;
	      q = tmp;
	    }	   
	    	    
	  while(!r->xBips)
	    {
	      if(!r->xBips)
		newviewBipartitions(r);
	    }	   

	  for(i = 0; i < vectorLength; i++)
	    vector[i] = left[i] | right[i];	    	 
	}
      else
	{	    
	  while((!r->xBips) || (!q->xBips))
	    {
	      if(!q->xBips)
		newviewBipartitions(q);
	      if(!r->xBips)
		newviewBipartitions(r);
	    }	   	    	    	    	   

	  for(i = 0; i < vectorLength; i++)
	    vector[i] = left[i] | right[i];	 
	}
    }     
}


/** 
    @brief adds all bipartitions in the tree into a frequency hash 

    @param traln -- the tree from which bipartitions should be extracted   
    @param slot -- the slot (e.g., chain id)
 */ 
void BipartitionHash::addBipartitionsToHash(TreeAln &traln, nat slot)
{
  tree *tr = traln.getTr();

  assert(slot < numSlots); 


  // check, if every tree that is thrown into this hash has the same taxa at the same position  
  if(taxonNames.size() == 0)	
    {
      taxonNames.push_back(""); 
      for(int i = 1; i < tr->mxtips+ 1 ; ++i)
	taxonNames.push_back(string(tr->nameList[i])); 
    }
  else
    {
      for(int i = 1; i < tr->mxtips + 1; ++i)
	assert(strcmp(tr->nameList[i], taxonNames[i].c_str() ) == 0); 
    }


  // assign random hash numbers to tree
  for(int i = 1; i < tr->mxtips + 1 ; ++i)
    tr->nodep[i]->hash = randomHash[i]; 

  resetBitVectors();
  
  int cnt = 0; 
  extractBipartitions(traln, tr->nodep[1]->back, &cnt, slot);
  assert(cnt == tr->mxtips -3); 
}




double BipartitionHash::averageDeviationOfSplitFrequencies(double ignoreFreq)
{
  double asdsf = 0; 

  int *numSampled = (int*)exa_calloc(numSlots, sizeof(int)); 
  
  entry** bvs = (entry**)exa_malloc(sizeof(entry*) * h->entryCount); 
  nat ctr = 0; 
  
  for(nat i = 0; i < h->tableSize; ++i)
    {
      entry *e = h->table[i]; 
      
      while(e != NULL)
	{
	  for(nat j = 0; j < numSlots; ++j)
	    {
	      numSampled[j] += e->treeVector[j]; 
	    }
	  
	  bvs[ctr] = e;
	  ctr++; 
	  e = e->next; 
	}
    }
  assert(ctr == h->entryCount); 

  /* asserting */
  for(nat i = 0; i < numSlots; ++i)
    assert(numSampled[0] == numSampled[i]); 
  assert(numSampled[0]  %  (numTax -3 )  == 0) ; 

  double treesSampled = numSampled[0] / (numTax - 3  ); 
#ifdef DEBUG_ASDSF_PRINT_ALL_BIPS
  printf("trees sampled: %g\n", treesSampled)  ;
#endif

  int cntRelevant = 0; 
  for(nat i = 0; i < h->entryCount ; ++i)
    {
      entry *e = bvs[i]; 

      boolean relevant = FALSE; 
      for(nat j  = 0; j < numSlots; ++j)
	if( e->treeVector[j] >  (ignoreFreq * treesSampled) )
	  relevant = TRUE; 

	  if(relevant)
	    {
	      double 
		n = numSlots,
		sd = 0; 

	      double mu = 0; 
	      for(int j = 0; j < n; ++j)
		mu += (double)e->treeVector[j] / treesSampled; 
	      mu /= n ; 

	      for(int j = 0; j < n; ++j)
		sd += pow((( ((double)e->treeVector[j]) / treesSampled) - mu ) ,2); 
	      sd = sqrt(sd/n) ; 

#ifdef DEBUG_ASDSF_PRINT_ALL_BIPS
	      printBv(e->bitVector, numTaxa); 
	      printf("\t");
	      for(int  j = 0;  j < n ; j++)
		printf("%d,", e->treeVector[j]); 
	      printf("\tmu=%f\tsd=%f\n", mu, sd); 
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

  exa_free(numSampled); 

  return asdsf ; 
}
