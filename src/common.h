/**
   @file common.h
   @brief Defines, typdes and developmental switches needed by almost every file.  
*/ 

#ifndef _COMMON_H
#define _COMMON_H

#include "config.h"


#if HAVE_PLL != 0
#define exa_realloc rax_realloc
#define exa_malloc_aligned  rax_malloc_aligned
#define exa_free rax_free
#define exa_malloc rax_malloc
#define exa_calloc rax_calloc
#else 
#define exa_malloc_aligned malloc_aligned
#define exa_realloc realloc
#define exa_free free
#define exa_malloc malloc
#define exa_calloc calloc
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


/* #define _GO_TO_TREE_MOVE_INTEGARTION */
/* #define _GO_TO_INTEGRATION_MODE */
/* #define _EXPERIMENTAL_INTEGRATION_MODE */

#define _DISABLE_INIT_LNL_CHECK
/* #define INTEGRATION_GENERATION 3e4 */
#define INTEGRATION_GENERATION 1e6

/* #define UNSURE  */
/* #define EFFICIENT  */

/* #define PRINT_EVAL_CHOICE */

/* debugging */
/* #define DEBUG_SHOW_EACH_PROPOSAL */
/* #define DEBUG_LNL_VERIFY */
/* #define DEBUG_VERIFY_LNPR */

#endif
