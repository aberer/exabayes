/**
   @file axml.h
   @brief A wrapper to include all files from pll/examl.   
*/ 


#ifndef _AXML_H_
#define _AXML_H_

#include "config.h"
#include "common.h"


#ifdef HAVE_AVX
#define __AVX 
#endif


#if HAVE_PLL != 0
#else 
#include "mpi.h"
#endif


#ifdef __cplusplus
extern "C"
{
#endif

#if HAVE_PLL != 0
#include "pll/axml-pll.h"
#else 
#include "examl/axml-examl.h"
#endif


#ifdef __cplusplus
}
#endif

#endif
