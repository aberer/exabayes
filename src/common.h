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


#define NO_SEC_BL_MULTI
#define NOT_IMPLEMENTED  0
#define TODO 0

typedef  unsigned int nat ; 

#define TARGET_RATIO 0.25    ///  the golden acceptance ratio, we want to achieve
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
/* #define CONTROL_ESPR */
/* #define DEBUG_TREE_LENGTH */
/* #define DEBUG_PARS_SPR */

/* for debugging:  */

#define DEBUG_GUIDED_SPR 0 	/* dont comment out, set to 0 for deactivation  */
/* #define DEBUG_SHOW_TREE */
#define DEBUG_SHOW_EACH_PROPOSAL
/* #define DEBUG_CHECK_TREE_CONSISTENCY */
/* #define DEBUG_SHOW_TOPO_CHANGES */
/* #define VERIFY_LNL_SUPER_EXPENSIVE */
/* #define DEBUG_ASDSF_PRINT_ALL_BIPS */



#define DEBUG_LNL_VERIFY
/* BAD BAD BAD  */
/* #define DEBUG_VERIFY_LNPR	/\* verify the log prior probability  *\/ */


#endif
