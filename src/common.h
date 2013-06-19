/**
   @file common.h
   @brief Defines, typdes and developmental switches needed by almost every file.  
*/ 

#ifndef _COMMON_H
#define _COMMON_H


#define NO_SEC_BL_MULTI


#define NOT_IMPLEMENTED  0
/* #define PRINT printBothOpen */


#include "config.h"

/* #define PRINT printBothOpen */

#if HAVE_PLL != 0
#define exa_realloc rax_realloc  
#define exa_malloc_aligned  rax_malloc_aligned
#define exa_malloc rax_malloc
#define exa_free rax_free
#define exa_calloc rax_calloc
#else 

#include <stdlib.h>
#define exa_realloc realloc
#define exa_malloc_aligned malloc_aligned
#define exa_free free 
#define exa_malloc malloc
#define exa_calloc calloc
#endif

#define TODO 0


// TODO: we also need a header cleanup, just imported this from examl 
#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
  special bug fix, enforces denormalized numbers to be flushed to zero,
  without this program is a tiny bit faster though.
  #include <emmintrin.h> 
  #define MM_DAZ_MASK    0x0040
  #define MM_DAZ_ON    0x0040
  #define MM_DAZ_OFF    0x0000
*/
#endif


#define ALL_MODELS -1 

typedef  unsigned int nat ; 

#define NOT ! 


#define TARGET_RATIO 0.234    ///  the golden acceptance ratio, we want to achieve
#define ACCEPTED_LIKELIHOOD_EPS 1e-6
#define ACCEPTED_LNPR_EPS 1e-6 


/* ABOVE: stuff that is needed by everyone and can be defined
   repeatedly

   BELOW: some switches (e.g., debugging) and constants that are could
   be user variables later (for development lets just keep them here).

 */

/* #define DEBUG_ACCEPTANCE */

/* #define DEBUG_HASTINGS */

/* #define UNSURE  */
/* #define EFFICIENT  */

/* #define INCORRECT */


/* #define PRINT_MULT */

/* #define DEBUG_EVAL */
/* #define DEBUG_ARRAY_SWAP */



/* #define  TUNE_ONLY_IF_ENOUGH	/// only tune a parameter, once TUNE_FREQUENCY times the respective function had been called     */
/* #define DEBUG_PRINT_TUNE_INFO */

/* #define ENABLE_PRSF */
#define WEIGHT_EPS 1e-3  	/* guided spr    */

/* #define CONTROL_ESPR */


/* #define DEBUG_TREE_LENGTH */

/* #define DEBUG_PARS_SPR */

/* for debugging:  */
#define STRETCH_FACTOR 2 
#define DEBUG_GUIDED_SPR 0 	/* dont comment out, set to 0 for deactivation  */
/* #define DEBUG_SHOW_TREE */
/* #define DEBUG_SHOW_EACH_PROPOSAL */
/* #define DEBUG_CHECK_TREE_CONSISTENCY */
/* #define DEBUG_SHOW_TOPO_CHANGES */
/* #define VERIFY_LNL_SUPER_EXPENSIVE */
/* #define DEBUG_ASDSF_PRINT_ALL_BIPS */

/* #define DEBUG_VERIFY_LNPR	/\* verify the log prior probability  *\/ */
/* #define DEBUG_LNL_VERIFY */

#endif
