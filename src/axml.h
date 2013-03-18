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

#if HAVE_PLL == 1 
#include "pll/axml-pll.h"
#else 
#include "examl/axml-examl.h"
#endif


#endif
