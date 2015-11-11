#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"
#include "phylip.h"
#include "ssort.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define SWAP(x,y) do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)

#define CONSUME(x)         while (token.class & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

//struct rawdata
// {
//   unsigned char ** oa;         /* original alignment */
//   unsigned char ** pats;       /* unique site patterns */
// };


static char * 
readFile (const char * filename, int * n)
{
  FILE * fp;
  char * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1) return (NULL);
  *n = ftell (fp);
  if (*n == -1) return (NULL);
  rewind (fp);

  rawdata = (char *) malloc ((*n) * sizeof (char));
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

  fclose (fp);

  return (rawdata);
}

static int
printTokens (int input)
{
  struct ltoken_t token;

  do
   {
     NEXT_TOKEN

     /* begin of parser */
     switch (token.class)
      {
        case LEX_NUMBER:
          printf ("LEX_NUMBER (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case LEX_STRING:
          printf ("LEX_STRING (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case LEX_EOF:
          printf ("LEX_EOF\n");
          break;
        case LEX_WHITESPACE:
          printf ("LEX_WHITESPACE\n");
          break;
        case LEX_NEWLINE:
          printf ("LEX_NEWLINE\n");
          break;
        case LEX_UNKNOWN:
          printf ("LEX_UNKNOWN (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
      }
     /* end of parser */


   }
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN);

  if (token.class == LEX_UNKNOWN) return (0);

  return (1);
}

static struct phylip_t *
alloc_phylip_struct (int nTaxa, int seqLen)
 {
   int i;
   struct phylip_t * phylip;
   void * mem;
   
   /** TODO */
   phylip = (struct phylip_t *) malloc (sizeof (struct phylip_t));
   phylip->seq = (unsigned char **) malloc ((nTaxa + 1) * sizeof (unsigned char *));
   mem = malloc (sizeof (unsigned char) * (seqLen + 1) * nTaxa);
   for (i = 1; i <= nTaxa; ++i)
    {
      phylip->seq[i] = (unsigned char *) (mem + (i - 1) * (seqLen + 1) * sizeof (unsigned char));
      phylip->seq[i][seqLen] = 0;
    }
   phylip->seq[0] = NULL;
    
   phylip->label = (char **) malloc ((nTaxa + 1) * sizeof (char *));

   phylip->nTaxa   = nTaxa;
   phylip->seqLen  = seqLen;
   phylip->weights = NULL;

   return (phylip);
 }

void
free_phylip_struct (struct phylip_t * phylip)
{
  int i;

  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     free (phylip->label[i]);
   }
  free (phylip->label);
  free (phylip->seq[1]);
  free (phylip->seq);
  free (phylip);
}


void 
dump_struct (struct phylip_t * pd)
 {
   int i;

   printf ("=> Dumping phylip_t\n");
   printf ("%d %d\n", pd->nTaxa, pd->seqLen);
   for (i = 0; i < pd->nTaxa; ++ i)
    {
      printf ("|%s| |%s|\n", pd->label[i], pd->seq[i]);
    }
 }

int
read_phylip_header (char * rawdata, int * inp, int * nTaxa, int * seqLen)
{
  struct ltoken_t token;
  int input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

  if (token.class != LEX_NUMBER) return (0);

  *nTaxa = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  if (token.class != LEX_NUMBER) return (0);

  *seqLen = atoi (token.lexeme);

  *inp = input;

  return (*nTaxa && *seqLen);
}

static inline int
parsedOk (int * actLen, int nTaxa, int seqLen  )
{
  int i;

  for (i = 1; i <= nTaxa; ++ i)
   {
     if (actLen[i] != seqLen) return (0);
   }
  
  return (1);
}


static int
parse_phylip (char * rawdata, struct phylip_t * phylip, int input)
{
  int i,j;
  struct ltoken_t token;
  int * seqLen;
  int rc;

  seqLen = (int *) calloc (phylip->nTaxa + 1, sizeof (int));

  NEXT_TOKEN
  for (i = 0; ; ++i)
  {
    j = i % phylip->nTaxa;
    if (i < phylip->nTaxa) 
     {
       if (token.class == LEX_EOF)
        {
          rc = parsedOk (seqLen, phylip->nTaxa, phylip->seqLen);
          free (seqLen);
          return (rc);
        }

       if (token.class == LEX_UNKNOWN)
        {
          free (seqLen);
          return (0);
        }

       CONSUME(LEX_WHITESPACE | LEX_NEWLINE)


       if (token.class != LEX_STRING && token.class != LEX_NUMBER)
        {
          free (seqLen);
          return (0);
        }
       phylip->label[i + 1] = strndup (token.lexeme, token.len);
       NEXT_TOKEN
       CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
     }
    
    while (1)
     {
       if (token.class == LEX_EOF)
        {
          rc = parsedOk (seqLen, phylip->nTaxa, phylip->seqLen);
          free (seqLen);
          return (rc);
        }

       if (token.class == LEX_UNKNOWN)
        {
         free (seqLen);
         return (0);
        }
       
       if (token.class == LEX_NEWLINE) break;

       if (token.class != LEX_STRING)
        {
          free (seqLen);
          return (0);
        }

       if (seqLen[j + 1] + token.len > phylip->seqLen) 
        {
          fprintf (stderr, "Sequence %d is larger than specified\n", j + 1);
          free (seqLen);
          return (0);
        }
       memmove (phylip->seq[j + 1] + seqLen[j + 1], token.lexeme, token.len);
       seqLen[j + 1] += token.len;

       NEXT_TOKEN
       CONSUME (LEX_WHITESPACE)
     }
    CONSUME(LEX_WHITESPACE | LEX_NEWLINE);
  }
}

/* Phylip parsers. Use the following attributed grammar 
 * 
 *        S -> HEADER ENDL DATA
 *   HEADER -> LEX_NUMBER LEX_WHITESPACE LEX_NUMBER ENDL |
 *             LEX_WHITESPACE LEX_NUMBER LEX_WHITESPACE LEX_NUMBER ENDL
 *     ENDL -> LEX_WHITESPACE LEX_NEWLINE | LEX_NEWLINE
 *     DATA -> LEX_STRING LEX_WHITESPACE LEX_STRING ENDL DATA |
 *             LEX_WHITESPACE LEX_STRING LEX_WHITESPACE LEX_STRING ENDL DATA | 
 *             LEX_STRING LEX_WHITESPACE LEX_STRING LEX_EOF |
 *             LEX_WHITESPACE LEX_STRING LEX_WHITESPACE LEX_STRING LEX_EOF
 */
struct phylip_t *
pllPhylipParse (const char * filename)
{
  int n, input, nTaxa, seqLen;
  char * rawdata;
  struct phylip_t * phylip;

  rawdata = readFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }

  init_lexan (rawdata, n);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (rawdata, &input, &nTaxa, &seqLen))
   {
     free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     return (0);
   }

  /* allocate the phylip structure */
  phylip = alloc_phylip_struct (nTaxa, seqLen);

  if (! parse_phylip (rawdata, phylip, input))
   {
     return (0);
     printf ("Finished with error in parsing ...\n");
     free_phylip_struct (phylip);
     free (rawdata);
     return (0);
   }
  
  free (rawdata);

  return (phylip);
}

void
pllPhylipExec (struct phylip_t * phylip, int op)
{
  void * mem;
  char ** sites;
  int i, j, k;
  int * oi;
  int dups = 0;


  sites = (char **) malloc (phylip->seqLen * sizeof (char *));
  
  mem = malloc ((phylip->nTaxa + 1) * phylip->seqLen * sizeof (char));
  for (i = 0; i < phylip->seqLen; ++ i)
   {
     sites[i] = (char *) (mem + i * (phylip->nTaxa + 1) * sizeof (char));
   }

  for (i = 0; i < phylip->seqLen; ++ i)
   {
     for (j = 0; j < phylip->nTaxa; ++ j)
      {
        sites[i][j] = phylip->seq[j + 1][i];
      }
     sites[i][j] = 0;
   }

  printf ("Original sequences:\n");
  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     printf ("%s\n", phylip->seq[i]);
   }

/*
  printf ("\nTransposed sequences\n");
  for (i = 0; i < phylip->seqLen; ++ i)
   {
     printf ("%s\n", sites[i]);
   }
*/
  oi = ssort1main (sites, phylip->seqLen);

/*
  printf ("\nSorted Transposed sequences\n");
  for (i = 0; i < phylip->seqLen; ++ i)
   {
     printf ("%s\n", sites[i]);
   }
*/

  for (i = 0; i < phylip->seqLen; ++ i) oi[i] = 1;

  for (i = 1; i < phylip->seqLen; ++ i)
   {
     if (! strcmp (sites[i], sites[i - 1]))
      {
        ++dups;
        oi[i] = 0;
      }
   }

  free (phylip->seq[1]);

  phylip->seqLen = phylip->seqLen - dups;
  mem = malloc ((phylip->seqLen + 1) * sizeof (unsigned char) * phylip->nTaxa);
  for (i = 0; i < phylip->nTaxa; ++ i)
   {
     phylip->seq[i + 1] = (unsigned char *) (mem + i * (phylip->seqLen + 1) * sizeof (unsigned char));
     phylip->seq[i + 1][phylip->seqLen] = 0;
   }

  for (i = 0, k = 0; i < phylip->seqLen + dups; ++ i)
   {
     if (!oi[i]) 
      {
        ++k;
        continue;
      }

     for (j = 0; j < phylip->nTaxa; ++ j)
      {
        phylip->seq[j + 1][i - k] = sites[i][j];
      }
   }

  printf ("\nUnique sequences:\n");
  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     printf ("%s\n", phylip->seq[i]);
   }


  printf ("\nTaxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  free (oi);

  free (sites);
}

//void pl_phylip_subst (struct phylip_t * pd, int type)

void
pllPhylipDump (struct phylip_t * phylip)
{
  FILE * fp;
  int i;

  fp = fopen ("output.seq.phy","w");
  fprintf (fp, "%d %d\n", phylip->nTaxa, phylip->seqLen);
  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     fprintf (fp, "%s %s\n", phylip->label[i], phylip->seq[i]);
   }

  fclose (fp);
}

void
pllPhylipSubst (struct phylip_t * phylip, int type)
{
  unsigned char meaningDNA[256];
  unsigned char  meaningAA[256];
  unsigned char * d;
  int i, j;

  for (i = 0; i < 256; ++ i)
   {
     meaningDNA[i] = -1;
     meaningAA[i]  = -1;
   }

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['Y'] = 10;
  meaningDNA['a'] =  1;
  meaningDNA['b'] = 14;
  meaningDNA['c'] =  2;
  meaningDNA['d'] = 13;
  meaningDNA['g'] =  4;
  meaningDNA['h'] = 11;
  meaningDNA['k'] = 12;
  meaningDNA['m'] =  3;
  meaningDNA['r'] =  5;
  meaningDNA['s'] =  6;
  meaningDNA['t'] =  8;
  meaningDNA['u'] =  8;
  meaningDNA['v'] =  7;
  meaningDNA['w'] =  9;
  meaningDNA['y'] = 10;

  meaningDNA['N'] =
  meaningDNA['n'] =
  meaningDNA['O'] =
  meaningDNA['o'] =
  meaningDNA['X'] =
  meaningDNA['x'] =
  meaningDNA['-'] =
  meaningDNA['?'] = 15;
 
  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/
  meaningAA['a'] =  0;  /* alanine */
  meaningAA['r'] =  1;  /* arginine */
  meaningAA['n'] =  2;  /*  asparagine*/
  meaningAA['d'] =  3;  /* aspartic */
  meaningAA['c'] =  4;  /* cysteine */
  meaningAA['q'] =  5;  /* glutamine */
  meaningAA['e'] =  6;  /* glutamic */
  meaningAA['g'] =  7;  /* glycine */
  meaningAA['h'] =  8;  /* histidine */
  meaningAA['i'] =  9;  /* isoleucine */
  meaningAA['l'] =  10; /* leucine */
  meaningAA['k'] =  11; /* lysine */
  meaningAA['m'] =  12; /* methionine */
  meaningAA['f'] =  13; /* phenylalanine */
  meaningAA['p'] =  14; /* proline */
  meaningAA['s'] =  15; /* serine */
  meaningAA['t'] =  16; /* threonine */
  meaningAA['w'] =  17; /* tryptophan */
  meaningAA['y'] =  18; /* tyrosine */
  meaningAA['v'] =  19; /* valine */
  meaningAA['b'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
  meaningAA['x'] = 
  meaningAA['?'] = 
  meaningAA['*'] = 
  meaningAA['-'] = 22;

  d = (type == PHYLIP_DNA_DATA) ? meaningDNA : meaningAA; 

  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     for (j = 0; j < phylip->seqLen; ++ j)
      {
        phylip->seq[i][j] = d[phylip->seq[i][j]];
      }
   }
}

void 
usage (const char * cmd_name)
{
  fprintf (stderr, "Usage: %s [PHYLIP-FILE]\n", cmd_name);
}

