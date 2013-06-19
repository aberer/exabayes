#include <cassert>

#include "axml.h"
#include "GlobalVariables.hpp"
#include "TreeAln.hpp"

static void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL != 0
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}

static void  treeEchoContext (FILE *fp1, FILE *fp2, int n)
{ /* treeEchoContext */
  int      ch;
  boolean  waswhite;
  
  waswhite = TRUE;
  
  while (n > 0 && ((ch = getc(fp1)) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = TRUE;
    }
    else {
      waswhite = FALSE;
    }
    
    if (ch > '\0') {putc(ch, fp2); n--;}
  }
}



static boolean treeLabelEnd (int ch)
{
  switch (ch) 
    {
    case EOF:  
    case '\0':  
    case '\t':  
    case '\n':  
    case '\r': 
    case ' ':
    case ':':  
    case ',':   
    case '(':   
    case ')':  
    case ';':
      return TRUE;
    default:
      break;
    }
  return FALSE;
}


static boolean  treeGetLabel (FILE *fp, char *lblPtr, int maxlen)
{
  int      ch;
  boolean  done, quoted, lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *) NULL; 
  else 
    if (lblPtr == NULL) 
      maxlen = 0;

  ch = getc(fp);
  done = treeLabelEnd(ch);

  lblfound = ! done;
  quoted = (ch == '\'');
  if (quoted && ! done) 
    {
      ch = getc(fp); 
      done = (ch == EOF);
    }

  while (! done) 
    {
      if (quoted) 
	{
	  if (ch == '\'') 
	    {
	      ch = getc(fp); 
	      if (ch != '\'') 
		break;
	    }
        }
      else 
	if (treeLabelEnd(ch)) break;     

      if (--maxlen >= 0) *lblPtr++ = ch;
      ch = getc(fp);
      if (ch == EOF) break;
    }

  if (ch != EOF)  (void) ungetc(ch, fp);

  if (lblPtr != NULL) *lblPtr = '\0';

  return lblfound;
}


static boolean  treeFlushLabel (FILE *fp)
{ 
  return  treeGetLabel(fp, (char *) NULL, (int) 0);
}





static int treeFinishCom (FILE *fp, char **strp)
{
  int  ch;
  
  while ((ch = getc(fp)) != EOF && ch != ']') {
    if (strp != NULL) *(*strp)++ = ch;    /* save character  */
    if (ch == '[') {                      /* nested comment; find its end */
      if ((ch = treeFinishCom(fp, strp)) == EOF)  break;
      if (strp != NULL) *(*strp)++ = ch;  /* save closing ]  */
    }
  }
  
  if (strp != NULL) **strp = '\0';        /* terminate string  */
  return  ch;
} /* treeFinishCom */



static int treeGetCh (FILE *fp)         /* get next nonblank, noncomment character */
{ /* treeGetCh */
  int  ch;

  while ((ch = getc(fp)) != EOF) {
    if (whitechar(ch)) ;
    else if (ch == '[') {                   /* comment; find its end */
      if ((ch = treeFinishCom(fp, (char **) NULL)) == EOF)  break;
    }
    else  break;
  }
  
  return  ch;
}



static boolean treeNeedCh (FILE *fp, int c1, string where)
{
  int  c2;
  
  if ((c2 = treeGetCh(fp)) == c1)  return TRUE;
  
  printf("ERROR: Expecting '%c' %s tree; found:", c1, where.c_str());
  if (c2 == EOF) 
    {
      printf("End-of-File");
    }
  else 
    {      	
      ungetc(c2, fp);
      treeEchoContext(fp, stdout, 40);
    }
  putchar('\n');

  if(c1 == ':')    
    printf("RAxML may be expecting to read a tree that contains branch lengths\n");

  return FALSE;
}






static int treeFindTipByLabelString(char  *str, tree *tr)                    
{
  int lookup = lookupWord(str, tr->nameHash);

  if(lookup > 0)
    {
      assert(! tr->nodep[lookup]->back);
      return lookup;
    }
  else
    { 
      printf("ERROR: Cannot find tree species: %s\n", str);
      return  0;
    }
}


static int treeFindTipName(FILE *fp, tree *tr)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabel(fp, str, nmlngth+2))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   

  return  n;
}


static boolean treeProcessLength (FILE *fp, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(fp)) == EOF)  return FALSE;    /*  Skip comments */
  (void) ungetc(ch, fp);
  
  if (fscanf(fp, "%lf", dptr) != 1) {
    printf("ERROR: treeProcessLength: Problem reading branch length\n");
    treeEchoContext(fp, stdout, 40);
    printf("\n");
    return  FALSE;
  }
  
  return  TRUE;
}



static int treeFlushLen (FILE  *fp)
{
  double  dummy;  
  int     ch;
  
  ch = treeGetCh(fp);
  
  if (ch == ':') 
    {
      ch = treeGetCh(fp);
      
      ungetc(ch, fp);
      if(!treeProcessLength(fp, & dummy)) return 0;
      return 1;	  
    }
  
  
  
  if (ch != EOF) (void) ungetc(ch, fp);
  return 1;
}


static boolean addElementLen (FILE *fp, tree *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount)
{   
  nodeptr  q;
  int      n, ch, fres;
  
  if ((ch = treeGetCh(fp)) == '(') 
    { 
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (tr->rooted || n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      return FALSE;
	    }
	  else 
	    {
	      assert(!readNodeLabels);
	      tr->rooted = TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (! addElementLen(fp, tr, q->next, readBranchLengths, readNodeLabels, lcount))        return FALSE;
      if (! treeNeedCh(fp, ',', "in"))             return FALSE;
      if (! addElementLen(fp, tr, q->next->next, readBranchLengths, readNodeLabels, lcount))  return FALSE;
      if (! treeNeedCh(fp, ')', "in"))             return FALSE;
      
      if(readNodeLabels)
	{
	  char label[64];
	  int support;

	  if(treeGetLabel (fp, label, 10))
	    {	
	      int val = sscanf(label, "%d", &support);
      
	      assert(val == 1);

	      /*printf("LABEL %s Number %d\n", label, support);*/
	      /*p->support = q->support = support;*/
	      /*printf("%d %d %d %d\n", p->support, q->support, p->number, q->number);*/
	      assert(p->number > tr->mxtips && q->number > tr->mxtips);
	      *lcount = *lcount + 1;
	    }
	}
      else	
	(void) treeFlushLabel(fp);
    }
  else 
    {   
      ungetc(ch, fp);
      if ((n = treeFindTipName(fp, tr)) <= 0)          return FALSE;
      q = tr->nodep[n];
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
    }
  
  if(readBranchLengths)
    {
      double branch;
      if (! treeNeedCh(fp, ':', "in"))                 return FALSE;
      if (! treeProcessLength(fp, &branch))            return FALSE;

      // TODO only works for 1 trbl
      hookup(p, q, &branch, 1 );
    }
  else
    {
      fres = treeFlushLen(fp);
      if(!fres) return FALSE;
      
      exa_hookupDefault(tr, p, q);
    }
  return TRUE;          
}




void myTreeReadLen(FILE *fp, tree *tr, boolean hasBL)
{
  boolean 
    readBranches = FALSE,  
    readNodeLabels = FALSE, 
    topologyOnly = FALSE; 

  if(hasBL)
    {
      readBranches = TRUE; 
    }
  else 
    {
      topologyOnly = TRUE; 
    }
    
  nodeptr  
    p;
  
  int      
    i, 
    ch, 
    lcount = 0; 

  for (i = 1; i <= tr->mxtips; i++) 
    {
      tr->nodep[i]->back = (node *) NULL; 
      /*if(topologyOnly)
	tr->nodep[i]->support = -1;*/
    }

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;

      /*if(topologyOnly)
	{
	  tr->nodep[i]->support = -2;
	  tr->nodep[i]->next->support = -2;
	  tr->nodep[i]->next->next->support = -2;
	  }*/
    }

  if(topologyOnly)
    tr->start       = tr->nodep[tr->mxtips];
  else
    tr->start       = tr->nodep[1];

  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;      
 

  // assert numbranches == 1  
  for(i = 0; i < 1; i++)
    tr->partitionSmoothed[i] = FALSE;
  
  tr->rooted      = FALSE;     

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');


  /* if(!topologyOnly) */
  /*   assert(readBranches == FALSE && readNodeLabels == FALSE); */
  
       
  if (! addElementLen(fp, tr, p, readBranches, readNodeLabels, &lcount))                 
    assert(0);
  if (! treeNeedCh(fp, ',', "in"))                
    assert(0);
  if (! addElementLen(fp, tr, p->next, readBranches, readNodeLabels, &lcount))
    assert(0);
  if (! tr->rooted) 
    {
      if ((ch = treeGetCh(fp)) == ',') 
	{ 
	  if (! addElementLen(fp, tr, p->next->next, readBranches, readNodeLabels, &lcount))
	    assert(0);	    
	}
      else 
	{                                    /*  A rooted format */
	  tr->rooted = TRUE;
	  if (ch != EOF)  (void) ungetc(ch, fp);
	}	
    }
  else 
    {      
      p->next->next->back = (nodeptr) NULL;
    }
  if (! treeNeedCh(fp, ')', "in"))                
    assert(0);

  if(topologyOnly)
    assert(!(tr->rooted && readNodeLabels));

  (void) treeFlushLabel(fp);
  
  if (! treeFlushLen(fp))                         
    assert(0);
 
  if (! treeNeedCh(fp, ';', "at end of"))       
    assert(0);

  tr->start = findAnyTip(p, tr->mxtips);    
 
  assert(tr->ntips == tr->mxtips);
}



/* I have no clue, why getc does not work here...
   this is not extensively tested
 */
boolean nextTreeHasBranchLength(FILE *fp)
{
  boolean foundBL = FALSE; 
  long oldPos = ftell(fp); 

  int c = 0; 
  while(( c = getc(fp)) != EOF)
    {
      if(c == ':')
	{
	  foundBL = TRUE; 
	  break; 
	}
      if(c == '\n')
	break; 
    }

  fseek(fp,oldPos, SEEK_SET ); 
  return foundBL;  
}




boolean readTreeWithOrWithoutBL(tree *tr, FILE *fh)
{  
  boolean hasBL = nextTreeHasBranchLength(fh); 

  myTreeReadLen(fh, tr, hasBL); 

  /* 
     TODO this is a horrible hack
     if the last char left after reading a tree is a newline, the next tree read would stop immediately. 

     Still: this only works, if we have one tree per line (something
     the user really should be able to accomplish).
   */  
  int c = getc(fh); 
  if(c != EOF && c != '\n')
    ungetc(c, fh); 
 
  return hasBL; 
}







// TODO HACK 
void traverseInitCorrect(nodeptr p, int *count, shared_ptr<TreeAln> traln )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;

  for( i = 0; i < traln->getNumBranches(); i++)
    {
      double val = traln->getBranchLength(p, i); 
      {
	double tmp = exp( - val  / tr->fracchange); 
	traln->setBranchLengthBounded(tmp , i,p); 
      }
    }
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitCorrect(q->back, count, traln);
	  q = q->next;
	} 
    }
}
