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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with 
 *  thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>


#include <mpi.h>

extern MPI_Comm comm; 


typedef unsigned int nat; 

#include "axml.h"

#ifndef __SIM_SSE3
#error "error: could not find sse3. This is either a programming error or your machine does not support sse3."
#endif



#if (defined(__SIM_SSE3) && !defined(__AVX))

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#define INTS_PER_VECTOR 4
#define LONG_INTS_PER_VECTOR 2
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128

#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>
#include <pmmintrin.h>

#define INTS_PER_VECTOR 8
#define LONG_INTS_PER_VECTOR 4
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd

#endif

static inline int states2pos (int states)
{
  switch(states)
    {
    case 4 :
      return 0; 
    case 20: 
      return 1; 
    default : assert(0); 
    }
}


extern const nat mask32[32]; 
/* vector-specific stuff */


extern double masterTime;

/************************************************ pop count stuff ***********************************************/


nat bitcount_32_bit(nat i)
{
  return ((nat) __builtin_popcount(i));
}

/* bit count for 64 bit integers */

inline nat bitcount_64_bit(unsigned long i)
{
  return ((nat) __builtin_popcountl(i));
}


/* bit count for 128 bit SSE3 and 256 bit AVX registers */

#if (defined(__SIM_SSE3) || defined(__AVX))
static inline unsigned int vectorPopcount(INT_TYPE v)
{
  unsigned long
    counts[LONG_INTS_PER_VECTOR] __attribute__ ((aligned (BYTE_ALIGNMENT)));

  int    
    i,
    sum = 0;
  
  VECTOR_STORE((CAST)counts, v);

  for(i = 0; i < LONG_INTS_PER_VECTOR; i++)
    sum += __builtin_popcountl(counts[i]);
 	     
  return ((unsigned int)sum);
}
#endif






/********************************DNA FUNCTIONS *****************************************************************/


static int checkerPars(tree *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->mxtips))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
	return group;

      group = checkerPars(tr, p->next->back);
      if(group != -9) 
	return group;

      group = checkerPars(tr, p->next->next->back);
      if(group != -9) 
	return group;

      return -9;
    }
}

static boolean tipHomogeneityCheckerPars(tree *tr, nodeptr p, int grouping)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(tr->constraintVector[p->number] != grouping) 
	return FALSE;
      else 
	return TRUE;
    }
  else
    {   
      return  (tipHomogeneityCheckerPars(tr, p->next->back, grouping) && tipHomogeneityCheckerPars(tr, p->next->next->back,grouping));      
    }
}

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->xPars || (s = s->next)->xPars)
    {
      p->xPars = s->xPars;
      s->xPars = 0;
    }

  assert(p->next->xPars || p->next->next->xPars || p->xPars);

}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, boolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->xPars)
    getxnodeLocal(p);  
  
  if(full)
    {
       if(q->number > maxTips) 
	 computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->xPars) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->xPars) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}







#if (defined(__SIM_SSE3) || defined(__AVX))

static void newviewParsimonyIterativeFast(tree *tr)
{    
  INT_TYPE
    allOne = SET_ALL_BITS_ONE;

  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 
  
  for(index = 4; index < count; index += 4)
    {       
      
      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  nat totalScore = 0 ;

	  int numByte = tr->partitionData[model].parsimonyLength * tr->partitionData[model].states * 2 * tr->mxtips; 
	  parsimonyNumber *ptrStart = tr->partitionData[model].parsVect,
	    *ptrEnd = ptrStart + numByte; 

	  unsigned long
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;

	  switch(states)
	    {
	    case 2:       
	      {
		/* check carefully, if we have trouble with the indices  */
		assert(0); 
		
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *this[2];

		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(nat i = 0; i < width; i += INTS_PER_VECTOR)
		  { 
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C,
		      v_A, v_C;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);		  		  		  		  
		    
		    v_N = VECTOR_BIT_OR(l_A, l_C);
		    
		    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += vectorPopcount(v_N);		  
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber *left[4];
		parsimonyNumber *right[4];
		parsimonyNumber *this[4];

		for(k = 0; k < 4; k++)
		  {

		    int numL = (width * 4 * qNumber) + width * k; 
		    int numR = (width * 4 * rNumber) + width * k; 
		    int numT = (width * 4 * pNumber) + width * k; 

		    left[k]  = ptrStart + numL;
		    right[k] = ptrStart + numR;
		    this[k]  = ptrStart + numT;

		  }

		for(nat iter = 0; iter < (nat)width; iter += INTS_PER_VECTOR)
		  {	 	  
		    /* printf("at %d / %d\t INTS=%d\n", iter, width, INTS_PER_VECTOR) ;  */

		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C, l_G, l_T,
		      v_A, v_C, v_G, v_T;	    	 

		    s_l = VECTOR_LOAD((CAST)(&left[0][iter]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][iter]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);

		    s_l = VECTOR_LOAD((CAST)(&left[1][iter]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][iter]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);

		    s_l = VECTOR_LOAD((CAST)(&left[2][iter]));
		    s_r = VECTOR_LOAD((CAST)(&right[2][iter]));
		    l_G = VECTOR_BIT_AND(s_l, s_r);
		    v_G = VECTOR_BIT_OR(s_l, s_r);

		    s_l = VECTOR_LOAD((CAST)(&left[3][iter]));
		    s_r = VECTOR_LOAD((CAST)(&right[3][iter]));
		    l_T = VECTOR_BIT_AND(s_l, s_r);
		    v_T = VECTOR_BIT_OR(s_l, s_r);

		    v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));	  	 	    	  
		    
		    VECTOR_STORE((CAST)(&this[0][iter]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&this[1][iter]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
		    VECTOR_STORE((CAST)(&this[2][iter]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
		    VECTOR_STORE((CAST)(&this[3][iter]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T))); 
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += vectorPopcount(v_N);	
		  }
	      }
	      break;
	    case 20:
	      {
		/* check carefully, if we have trouble with the indices  */
		/* assert(0);  */

		parsimonyNumber
		  *left[20],
		  *right[20],
		  *this[20];

		for(k = 0; k < 20; k++)
		  {		    
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(nat i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[20], 
		      v_A[20];	    	 
		    
		    for(j = 0; j < 20; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < 20; j++)		    
		      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += vectorPopcount(v_N);
		  }
	      }
	      break;
	    default:
	      {
		/* check carefully, if we have trouble with the indices  */
		assert(0); 

		parsimonyNumber
		  *left[32], 
		  *right[32],
		  *this[32];

		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }

		for(nat i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[32], 
		      v_A[32];	    	 
		    
		    for(j = 0; j < states; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < states; j++)		    
		      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += vectorPopcount(v_N);
		  }	  			
	      }
	    }

	  tr->parsimonyScore[pNumber * tr->NumberOfModels + model ]
	    = totalScore + tr->parsimonyScore[rNumber * tr->NumberOfModels + model] 
	    + tr->parsimonyScore[qNumber * tr->NumberOfModels + model]; 	      
	}
    }
}

static void evaluateParsimonyIterativeFast(tree *tr,  nat *state2parsimony)
{
  INT_TYPE 
    allOne = SET_ALL_BITS_ONE;

  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;
 
  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr);
  
  for(nat i = 0; i < tr->NumberOfModels; ++i)
    {
      nat states = tr->partitionData[i].states; 
      state2parsimony[ states2pos( states ) ] = tr->parsimonyScore[pNumber * tr->NumberOfModels + i] + tr->parsimonyScore[qNumber * tr->NumberOfModels + i]; 
    }

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      nat
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     

	     nat sum = 0; 
	     
	     for(nat i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                       
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),		 
		   v_N = VECTOR_BIT_OR(l_A, l_C);
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += vectorPopcount(v_N); ;
	       }
	     /* pLengthAtBranch[model] = sum;  */
	     state2parsimony[states2pos(states)] += sum; 
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       *left[4],
	       *right[4];
      
	     parsimonyNumber *start = tr->partitionData[model].parsVect;

	     for(k = 0; k < 4; k++)
	       {
		 int numL = (width * 4 * qNumber) + width * k,
		   numR = (width * 4 * pNumber) + width * k; 
		 
		 left[k]  = start + numL;
		 right[k] = start + numR;
	       }        

	     nat sum = 0; 

	     for(nat i = 0; i < width; i += INTS_PER_VECTOR)
	       {

		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),
		   l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[2][i])), VECTOR_LOAD((CAST)(&right[2][i]))),
		   l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[3][i])), VECTOR_LOAD((CAST)(&right[3][i]))),
		   v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);

		 sum += vectorPopcount(v_N); 
	       }

	     /* pLengthAtBranch[model] = sum; */
	     state2parsimony[states2pos( states)] += sum; 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   

	      nat sum = 0; 

	      for(nat i = 0; i < width; i += INTS_PER_VECTOR)
		{                	       
		  int 
		    j;
		  
		  INT_TYPE      
		    l_A,
		    v_N = SET_ALL_BITS_ZERO;     
		  
		  for(j = 0; j < 20; j++)
		    {
		      l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		      v_N = VECTOR_BIT_OR(l_A, v_N);
		    }
		  
		  v_N = VECTOR_AND_NOT(v_N, allOne);
	
		  sum += vectorPopcount(v_N); ;
		}
	      /* pLengthAtBranch[model] = sum; */
	      state2parsimony[states2pos(states)] += sum; 
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       *left[32],  
	       *right[32]; 

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     nat sum = 0; 

	     for(nat i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	       
		 size_t
		   j;
		 
		 INT_TYPE      
		   l_A,
		   v_N = SET_ALL_BITS_ZERO;     
		 
		 for(j = 0; j < states; j++)
		   {
		     l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		     v_N = VECTOR_BIT_OR(l_A, v_N);
		   }
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);

		 sum += vectorPopcount(v_N); 
	       }
	     /* pLengthAtBranch[model] = sum;  */
	     state2parsimony[states2pos(states)] += sum; 
	   }
	 }
    }
}


#else

void newviewParsimonyIterativeFast(tree *tr)
{    
  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *this[2];
		
		parsimonyNumber
		   o_A,
		   o_C,
		   t_A,
		   t_C,	
		   t_N;
		
		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];		   

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		  
		    t_N = ~(t_A | t_C);	  

		    this[0][i] = t_A | (t_N & o_A);
		    this[1][i] = t_C | (t_N & o_C);		   
		    
		    totalScore += populationCount(t_N);   
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *this[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		parsimonyNumber
		   o_A,
		   o_C,
		   o_G,
		   o_T,
		   t_A,
		   t_C,
		   t_G,
		   t_T,	
		   t_N;

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];
		    t_G = left[2][i] & right[2][i];	  
		    t_T = left[3][i] & right[3][i];

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		    o_G = left[2][i] | right[2][i];	  
		    o_T = left[3][i] | right[3][i];

		    t_N = ~(t_A | t_C | t_G | t_T);	  

		    this[0][i] = t_A | (t_N & o_A);
		    this[1][i] = t_C | (t_N & o_C);
		    this[2][i] = t_G | (t_N & o_G);
		    this[3][i] = t_T | (t_N & o_T); 
		    
		    totalScore += populationCount(t_N);   
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *this[20];

		parsimonyNumber
		  o_A[20],
		  t_A[20],	  
		  t_N;

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    size_t k;
		    
		    t_N = 0;

		    for(k = 0; k < 20; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;

		    for(k = 0; k < 20; k++)		      
		      this[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += populationCount(t_N); 
		  }
	      }
	      break;
	    default:
	      {		
		parsimonyNumber
		  *left[32],
		  *right[32],
		  *this[32];
		
		parsimonyNumber
		  o_A[32],
		  t_A[32],	  
		  t_N;
		
		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }
		
		for(i = 0; i < width; i++)
		  {	 	  
		    t_N = 0;
		    
		    for(k = 0; k < states; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;
		    
		    for(k = 0; k < states; k++)		      
		      this[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += populationCount(t_N); 
		  }
	      }			      
	    } 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}



unsigned int evaluateParsimonyIterativeFast(tree *tr,  nat *state2parsimony)
{
  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber 
	       t_A,
	       t_C,	      
	       t_N,
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i++)
	       {                	                       
		 t_A = left[0][i] & right[0][i];
		 t_C = left[1][i] & right[1][i];
		 
		  t_N = ~(t_A | t_C);

		  sum += populationCount(t_N);    
		 
		 /* if(sum >= bestScore) */
		 /*   return sum;	 */	   	       
	       }
	     state2parsimony[states2pos(states)] += sum; 
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       t_A,
	       t_C,
	       t_G,
	       t_T,
	       t_N,
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i++)
	       {                	                        
		  t_A = left[0][i] & right[0][i];
		  t_C = left[1][i] & right[1][i];
		  t_G = left[2][i] & right[2][i];	  
		  t_T = left[3][i] & right[3][i];

		  t_N = ~(t_A | t_C | t_G | t_T);

		  sum += populationCount(t_N);     
		 
		 /* if(sum >= bestScore)		  */
		 /*   return sum;	         */
	       }	   	 
	     state2parsimony[states2pos(states)] += sum; 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i++)
		{ 
		  t_N = 0;
		  
		  for(k = 0; k < 20; k++)
		    {
		      t_A = left[k][i] & right[k][i];
		      t_N = t_N | t_A;
		    }
  	       
		  t_N = ~t_N;

		  sum += populationCount(t_N);      
		  
		  /* if(sum >= bestScore)	     */
		  /*   return sum;		    	        */
		}
	      state2parsimony[states2pos(states)] += sum; 
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[32], 
	       *right[32];  

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i++)
	       {                	       
		 t_N = 0;
		  
		 for(k = 0; k < states; k++)
		   {
		     t_A = left[k][i] & right[k][i];
		     t_N = t_N | t_A;
		   }
  	       
		  t_N = ~t_N;

		  sum += populationCount(t_N);      
		  		  		 
		 if(sum >= bestScore)			  
		   return sum;			   
	       }	     	     
	     state2parsimony[states2pos(states)] += sum; 
	   }
	 }
    }
  
  return sum;
}

#endif


void evaluateParsimony(tree *tr, nodeptr p, boolean full, nat *state2parsimony )
{
  /* volatile nat   */
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->xPars)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->xPars)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }

  ti[0] = counter;

  evaluateParsimonyIterativeFast(tr, state2parsimony);
}


void newviewParsimony(tree *tr,  nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(tr);
  }
}





/****************************************************************************************************************************************/


static boolean isInformative(tree *tr, int dataType, int site)
{
  assert(0);
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const nat
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);	           
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
	informativeCounter++;    
    } 
	  
  if(informativeCounter <= 1)
    return FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(check[j] > 1)
	    return TRUE;
	} 
    }
     
  return FALSE;	     
}


static void determineUninformativeSites(tree *tr,  int *informative)
{
  int 
    model,
    number = 0,
    i;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */


  for(model = 0; model < tr->NumberOfModels; model++)
    {
      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	{
	   if(isInformative(tr, tr->partitionData[model].dataType, i))
	     informative[i] = 1;
	   else
	     {
	       informative[i] = 0;
	       number++;
	     }  
	}      
    }

  /* printf("Uninformative Patterns: %d\n", number); */
}


static void compressDNA(tree *tr,  int *informative)
{
  size_t
    totalNodes,
    i,
    model;
   
  totalNodes = 2 * (size_t)tr->mxtips;

  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = (size_t)tr->partitionData[model].states,
	compressedEntries,
	compressedEntriesPadded,
	entries = 0, 
	lower = tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper;

      parsimonyNumber 
	**compressedTips = (parsimonyNumber **)exa_calloc(states,  sizeof(parsimonyNumber*)),
	*compressedValues = (parsimonyNumber *)exa_calloc(states,  sizeof(parsimonyNumber));
      
      for(i = lower; i < upper; i++)    
	if(informative[i])
	  entries += (size_t)tr->aliaswgt[i];     
  
      compressedEntries = entries / PCF;

      if(entries % PCF != 0)
	compressedEntries++;

#if (defined(__SIM_SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
	compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
	compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      
      int numByte = (size_t)compressedEntriesPadded * states * totalNodes; 

      tr->partitionData[model].parsVect = (parsimonyNumber *)exa_malloc_aligned(numByte * sizeof(parsimonyNumber));
     
      for(i = 0; i < numByte; i++)      
	tr->partitionData[model].parsVect[i] = 0;

      for(i = 0; i < (size_t)tr->mxtips; i++)
	{
	  size_t
	    w = 0,
	    compressedIndex = 0,
	    compressedCounter = 0,
	    index = 0;

	  for(k = 0; k < states; k++)
	    {
	      compressedTips[k] = &(tr->partitionData[model].parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
	      compressedValues[k] = 0;
	    }                
	      
	  for(index = lower; index < (size_t)upper; index++)
	    {
	      if(informative[index])
		{
		  const nat 
		    *bitValue = getBitVector(tr->partitionData[model].dataType);

		  parsimonyNumber 
		    value = bitValue[tr->yVector[i + 1][index]];	  
	      
		  for(w = 0; w < (size_t)tr->aliaswgt[index]; w++)
		    {	   
		      for(k = 0; k < states; k++)
			{
			  if(value & mask32[k])
			    compressedValues[k] |= mask32[compressedCounter];
			}
		     
		      compressedCounter++;
		  
		      if(compressedCounter == PCF)
			{
			  for(k = 0; k < states; k++)
			    {
			      compressedTips[k][compressedIndex] = compressedValues[k];
			      compressedValues[k] = 0;
			    }			 
			  
			  compressedCounter = 0;
			  compressedIndex++;
			}
		    }
		}
	    }
                           
	  for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	    {	
	      for(;compressedCounter < PCF; compressedCounter++)	      
		for(k = 0; k < states; k++)
		  compressedValues[k] |= mask32[compressedCounter];		  
	  
	      for(k = 0; k < states; k++)
		{
		  compressedTips[k][compressedIndex] = compressedValues[k];
		  compressedValues[k] = 0;
		}	      	      
	      
	      compressedCounter = 0;
	    }	 	
	}               
  
      tr->partitionData[model].parsimonyLength = compressedEntriesPadded;
      
      /* for test  */
      assert(0); 

      exa_free(compressedTips);
      exa_free(compressedValues);

      for(int i = 0; i < numByte; ++i)
	printf("%u,", tr->partitionData[model].parsVect[i]); 
      printf("\n");
    }

  tr->parsimonyScore = (nat*)exa_malloc_aligned(sizeof(nat) * totalNodes * tr->NumberOfModels);  
          
  for(i = 0; i < totalNodes * tr->NumberOfModels; i++) 
    tr->parsimonyScore[i] = 0;
}



static void stepwiseAddition(tree *tr,  nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;
  
  nat 
    mp;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
  tr->ti[1] = p->number;
  tr->ti[2] = p->back->number;

  nat state2parsimony[2] = {0,0}; 
  evaluateParsimonyIterativeFast(tr,  state2parsimony);

  mp = 0; 
  for(int i = 0; i < 2; ++i)
    mp += state2parsimony[i]; 

  MPI_Allreduce(MPI_IN_PLACE, &mp,1,MPI_UNSIGNED, MPI_SUM, comm);

  if(mp < tr->bestParsimony)
    {    
      tr->bestParsimony = mp;
      tr->insertNode = q;     
    }
 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips )
    {	      
      stepwiseAddition(tr,  p, q->next->back);
      stepwiseAddition(tr, p, q->next->next->back);
    }
}



void allocateParsimonyDataStructures(tree *tr)
{
  int 
    i,
    *informative = (int *)exa_malloc(sizeof(int) * (size_t)tr->originalCrunchedLength);
 
  determineUninformativeSites(tr, informative);

  compressDNA(tr, informative);

  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 1; i++)
    {
      nodeptr 
	p = tr->nodep[i];

      p->xPars = 1;		/*  */
      p->next->xPars = 0;
      p->next->next->xPars = 0;
    }

  tr->ti = (int*)exa_malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  exa_free(informative); 
}

void freeParsimonyDataStructures(tree *tr)
{
  size_t 
    model;

  exa_free(tr->parsimonyScore);
  
  for(model = 0; model < (size_t) tr->NumberOfModels; ++model)
    exa_free(tr->partitionData[model].parsVect);
  
  exa_free(tr->ti);
}


static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, 0.65);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
}

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, 0.65);
  hookupDefault(p->next->next, r, 0.65);
   
  newviewParsimony(tr, p);
}


static void buildSimpleTree (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], 0.65); /* whatever */
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, s, p);
}



static void makePermutationFast(int *perm, int n, tree *tr, nat seed)
{    
  tr->randomNumberSeed = seed; 
  int  
    i, 
    j, 
    k;

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {      
      double d =  randum(&tr->randomNumberSeed);

      k =  (int)((double)(n + 1 - i) * d);
      
      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

void makeParsimonyTreeFast(tree *tr, nat seed)
{   
  nodeptr  
    p = NULL, 
    f = NULL;    

  int 
    i, 
    nextsp,
    *perm        = (int *)exa_calloc((size_t)(tr->mxtips + 1) , sizeof(int));  

  unsigned int 
    randomMP, 
    startMP;         
  
  assert(!tr->constrained);

  makePermutationFast(perm, tr->mxtips, tr, seed);

  tr->ntips = 0;    
  
  tr->nextnode = tr->mxtips + 1;       
  
  buildSimpleTree(tr, perm[1], perm[2], perm[3]);

  f = tr->start;       
  
  while(tr->ntips < tr->mxtips) 
    {	
      nodeptr q;
      
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
        
      if(tr->grouped)
	{
	  int 
	    number = p->back->number;	  	 

	  tr->constraintVector[number] = -9;
	}
      
      stepwiseAddition(tr, q, f->back);

      /* printf("[%d / %d] bestnode: %d\n", tr->ntips, tr->mxtips, tr->insertNode->number);  */

      {
	nodeptr	  
	  r = tr->insertNode->back;
	
	int counter = 4;
	
	hookupDefault(q->next, tr->insertNode, 0.65);
	hookupDefault(q->next->next, r, 0.65);
	
	computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, FALSE);              
	tr->ti[0] = counter;
	
	newviewParsimonyIterativeFast(tr);
      }
    }    

  
  MPI_Barrier(comm);

  /* printf("ADD: %d\n", tr->bestParsimony);  */

  /* nodeRectifierPars(tr); */

  randomMP = tr->bestParsimony;        

  /* TODO removed the spr moves... does not work otherwise */

#if 0
  do
    {
      startMP = randomMP;
      nodeRectifierPars(tr);
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{
	  rearrangeParsimony(tr, pr, tr->nodep[i], 1, 20, FALSE);
	  
	  if(tr->bestParsimony < randomMP)
	    {		
	      restoreTreeRearrangeParsimony(tr, pr);
	      randomMP = tr->bestParsimony;
	    }
	}      		  	   
    }
  while(randomMP < startMP); 
#endif

  /* andre: i do not know, why the following line is necesary...   */
  /* nodeRectifierPars(tr); */
  /* tr->start = tr->nodep[1]; */
}
