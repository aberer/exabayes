#ifndef PARSER_H
#define PARSER_H
#include "phylip.h"


struct phylip_t * pllPhylipParse (const char *);
void pllPhylipExec (struct phylip_t **, int);


struct model_t * pllParseModel (const char * filename);


pllMerge (tree * tr, phylip_t * phylip, model_t * model);
#endif
