/**
   @file common.h
   @brief Defines, typdes and developmental switches needed by almost every file.  
*/ 

#ifndef _COMMON_H
#define _COMMON_H

#include "config.h"

/* only to be disabled for benchmarking! */
#define USE_NONBLOCKING_COMM

#define USE_SSE ( defined(HAVE_SSE3) &&  ! defined(MANUAL_SSE_OVERRIDE) ) 
#define USE_AVX ( defined(HAVE_AVX) && ! defined(MANUAL_AVX_OVERRIDE) && ! defined(MANUAL_SSE_OVERRIDE) )

#if USE_AVX
#define EXA_ALIGN VectAlign::AVX 
#elif USE_SSE
#define EXA_ALIGN VectAlign::SSE
#else 
#define EXA_ALIGN VectAlign::Normal
#endif

#define SOME_SCI_PRECISION std::scientific << std::setprecision(2)   
#define MAX_SCI_PRECISION  std::scientific << std::setprecision(std::numeric_limits<double>::digits10)   
#define SOME_FIXED_PRECISION std::fixed << std::setprecision(2)   
#define MORE_FIXED_PRECISION std::fixed << std::setprecision(4)   
#define PERC_PRECISION std::fixed << std::setprecision(2) 

#define NO_SEC_BL_MULTI
#define NOT_IMPLEMENTED  0
#define TODO 0

typedef unsigned int nat; 

#define TARGET_RATIO 0.25    ///  the golden acceptance ratio, we want to achieve
#define ACCEPTED_LIKELIHOOD_EPS 1e-6
#define ACCEPTED_LNPR_EPS 1e-6

/* some global switches */ 

/* #define _EXPERIMENTAL_INTEGRATION_MODE */

#define _DISABLE_INIT_LNL_CHECK

#define SHOW(sym) #sym << "=" << sym << "\t"

/* #define EFFICIENT  */

/* #define PRINT_EVAL_CHOICE */

/* debugging */
/* #define DEBUG_SHOW_EACH_PROPOSAL */
/* #define DEBUG_LNL_VERIFY */
/* #define DEBUG_VERIFY_LNPR */

/* many print statements  */
/* #define LNL_PRINT_DEBUG */

#endif
