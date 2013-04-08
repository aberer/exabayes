/**
   @file common.h
   @brief Defines, typdes and developmental switches needed by almost every file.  
*/ 

#ifndef _COMMON_H
#define _COMMON_H

#include "config.h"


#if HAVE_PLL == 1
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


typedef  unsigned int nat ; 

#define NOT ! 


#define INIT_BRANCH_LENGTHS 0.65
#define INIT_BL_MULT 1.386294
#define INIT_BL_SLID_WIN  0.075
#define INIT_RATE_SLID_WIN  0.15 
#define INIT_FREQ_SLID_WIN  0.2 
#define INIT_GAMMA_SLID_WIN  0.75
#define INIT_ESPR_MULT 0.098
#define INIT_GUIDED_RADIUS 5

#define TARGET_RATIO 0.234    ///  the golden acceptance ratio, we want to achieve
#define ACCEPTED_LIKELIHOOD_EPS 1e-6


/* ABOVE: stuff that is needed by everyone and can be defined
   repeatedly

   BELOW: some switches (e.g., debugging) and constants that are could
   be user variables later (for development lets just keep them here).

 */


#define ESPR_MULTIPLY_BL	 /// applyl multiplier to espr as by lakner

#define PRINT_FREQUENCY 500 
#define TUNE_FREQUENCY 100 

#define TUNE_PARAMETERS		/// turn off autotuning

/* #define DEBUG_PRINT_TUNE_INFO */


/* for debugging:  */
#define DEBUG_GUIDED_SPR 0 	/* dont comment out, set to 0 for deactivation  */
/* #define DEBUG_SHOW_TREE */
/* #define DEBUG_SHOW_EACH_PROPOSAL */
/* #define DEBUG_LNL_VERIFY */
/* #define DEBUG_CHECK_TOPO_CONSISTENCY 	/\* TODO <= merge more stuff into that *\/ */
/* #define DEBUG_SHOW_TOPO_CHANGES */


#define LENGTH_LNL_ARRAY 16 /// factor we need to multiply to the conditional arrays

/* #define DEBUG_ASDSF_PRINT_ALL_BIPS */

#define ASDSF_FREQ_CUTOFF 0.1	/*  ignore clades for which the frequency in no chain exceeds this value */
#define ASDSF_CONVERGENCE_CRITERION 0.005 ///  indicate convergence, as soon as the asdsf is below this value: 1-5% is considered good, 0.005 can be considered very good convergence 

/* a lot of debug information for the asdsf */


/* MC3 stuff for development */
#define HEAT_FACTOR 0.1
#define SWITCH_AFTER_GEN 1  	/* number of generations after  which we try to switch states   */
#define MC3_SPACE_FOR_TIME


/* here you can set the burnin 
   
   TODO make this a parameter in the config file  
 */
#define BURNIN 5000

#endif
