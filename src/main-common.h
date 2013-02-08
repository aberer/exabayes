#ifndef _MAIN_COMMON_H 
#define _MAIN_COMMON_H 

#include <string.h> 
#include <stdlib.h>
#include <stdio.h>




/* provides  */
int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg); 
void get_args(int argc, char *argv[], analdef *adef, tree *tr); 

#endif
