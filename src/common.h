#ifndef _COMMON_H
#define _COMMON_H

#include "config.h"


#if HAVE_PLL == 1
#define realloc rax_realloc  
#define malloc_aligned  rax_malloc_aligned
#define free rax_free
#define calloc rax_calloc
#else 
#include <stdlib.h>
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





/* ABOVE: stuff that is needed by everyone and can be defined
   repeatedly

   BELOW: some switches (e.g., debugging) and constants that are could
   be user variables later (for development lets just keep them here).

 */




#define _USE_NCL_PARSER
#define INIT_BRANCH_LENGTHS 0.65

/* for debugging:  */
/* #define DEBUG_SHOW_TREE */
/* #define DEBUG_SHOW_EACH_PROPOSAL */
// #define DEBUG_LNL_VERIFY



#define ASDSF_FREQ_CUTOFF 0.1	/*  ignore clades for which the frequency in no chain exceeds this value */
#define ASDSF_CONVERGENCE_CRITERION 0.05 /* indicate convergence, as soon as the asdsf is below this value  */

/* a lot of debug information for the asdsf */
#define ASDSF_BE_VERBOSE 

#endif
