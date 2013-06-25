/**
   @file axml.h
   @brief A wrapper to include all files from pll/examl.   
*/ 


#ifndef _AXML_H_
#define _AXML_H_

#include "config.h"
#include "common.h"

#ifdef _USE_GOOGLE_PROFILER
#include <google/profiler.h>
#endif





#ifdef HAVE_AVX
#define __AVX 
#endif


#if HAVE_PLL != 0

#else 



#ifdef _SEQUENTIAL

#ifdef __cplusplus
extern "C"{
#endif
#include "examl/mpiMock.h"
#ifdef __cplusplus
}
#endif
#else 
#include "mpi.h"
#endif

#endif


#ifdef __cplusplus
extern "C"{
#endif

#if HAVE_PLL != 0
#include "pll/axml-pll.h"
#else 
#include "examl/axml-examl.h"
#endif



#if HAVE_PLL == 0
  void newviewParsimony(tree *tr, nodeptr  p); 
  void evaluateParsimony(tree *tr, nodeptr p, boolean full, nat *partitionParsimony); 
#else 
  void newviewParsimony(tree *tr, partitionList *pr, nodeptr  p); 
  void evaluateParsimony(tree *tr, partitionList *pr, nodeptr p, boolean full, nat *partitionParsimony); 
#endif

#ifdef __cplusplus
}
#endif

#endif



