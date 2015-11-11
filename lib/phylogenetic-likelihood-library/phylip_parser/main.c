#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"
#include "phylip.h"
//#include "xalloc.h"
//#include "msa_sites.h"

int 
main (int argc, char * argv[])
{
  struct phylip_t * phylip;

  if (argc != 2)
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }
  
  phylip = pllPhylipParse (argv[1]);
  if (!phylip) 
   {
     printf ("Error while parsing\n");
     return (EXIT_FAILURE);
   }

  //dump_struct (pd);
  printf ("Taxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  
//  weight  = pl_phylip_deldups (&pd);
//  printf ("=> Eliminating dups\n");
//  dump_struct (pd);
//  for (i = 0; i < pd->seqlen; ++ i)
//  printf ("%d ", weight[i]);
//  printf ("\n");
  
  pllPhylipDump (phylip);
  //pllPhylipExec (phylip, 1);
  
  free_phylip_struct (phylip);
//  free (weight);

  /*
//  dump_struct (pd);

//  printf ("=> Sorting\n");
  // ms = construct_msa_sites (pd, SITES_CREATE | SITES_ELIMINATE_DUPLICATES | SITES_COMPUTE_WEIGHTS);
  ms = construct_msa_sites (pd, SITES_CREATE | SITES_COMPUTE_WEIGHTS);
  dump_sites (ms);


  for (i = 0; i < ms->seqlen; ++ i)
  printf ("%d ", ms->weight[i]);
  printf ("\n");

  free_phylip_struct (pd);
  pd = transpose (ms);
  free_sites_struct (ms);

  dump_struct (pd);
  for (i = 0; i < pd->seqlen; ++ i)
  printf ("%d ", pd->weight[i]);
  printf ("\n");
  free_phylip_struct (pd);
*/


  return (EXIT_SUCCESS);
}
