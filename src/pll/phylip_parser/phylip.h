#ifndef PHYLIP_H
#define PHYLIP_H

#define PHYLIP_KEEP_UNIQUE      0x00000001
#define PHYLIP_LEX_SORT         0x00000002
#define PHYLIP_SITE_WEIGHTS     0x00000004
#define PHYLIP_DNA_DATA         0x00000008
#define PHYLIP_PROT_DATA        0x00000010


struct phylip_t
 {
   int              nTaxa;
   int              seqLen;
   char          ** label;
   unsigned char ** seq;
   int            * weights;
 };

static struct phylip_t * alloc_phylip_struct (int taxa, int seqlen);
struct phylip_t * pllPhylipParse (const char *);
void pllPhylipExec (struct phylip_t *, int);
void pl_phylip_subst (struct phylip_t *, int);
void free_phylip_struct (struct phylip_t *);
void usage (const char * cmd_name);
void pllPhylipDump (struct phylip_t * phylip);

#endif
