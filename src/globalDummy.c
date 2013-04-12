

/**
   @brief this file contains legacy global variables. 
   
   It is compiled as a c-file such we have c-linkage here.
   
 */ 

#include "axml.h" 

#if HAVE_PLL == 1
#define GLOBAL_VARIABLES_DEFINITION
#endif

#define _INCLUDE_DEFINITIONS
#include "globalVariables.h"  
#undef  _INCLUDE_DEFINITIONS

#if HAVE_PLL == 1
#undef GLOBAL_VARIABLES_DEFINITION
#endif


/* inline unsigned int bitcount_64_bit(unsigned long i) */
/* { */
/*   printf("hello world!\n");  */
/* return 0 ;  */
/* } */


