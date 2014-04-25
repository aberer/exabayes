/**
   @file axml.h
   @brief A wrapper to include all files from pll/examl.   
*/ 


#ifndef _AXML_H_
#define _AXML_H_


#include "common.h"

/* extern int NUM_BRANCHES;  */

#ifdef _USE_GOOGLE_PROFILER
#include <google/profiler.h>
#endif



/* translate the defines for correct vectorization  */
#if USE_SSE
#define __SSE3
#endif

#if USE_AVX
#define __AVX 
#endif

#ifdef __cplusplus
extern "C"{
#endif

#include "lib/pll/pll-renamed.h"
#include "lib/pll/pllInternal.h"

unsigned int evaluateParsimony(pllInstance *tr, partitionList *pr, nodeptr p, boolean full );
void newviewParsimony(pllInstance *tr, partitionList *pr, nodeptr  p); 

#ifdef __cplusplus
}
#endif

#endif
