#include <cassert>
#include <sstream> 
#include <algorithm>

#include "axml.h"
#include "GlobalVariables.hpp"


class TreeAln; 

static void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL != 0
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}

static void  treeEchoContext ( std::istringstream &iss , int n)
{
  int      ch;
  boolean  waswhite;
  
  waswhite = TRUE;
  
  while (n > 0 && ((ch = iss.get()) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = TRUE;
    }
    else {
      waswhite = FALSE;
    }
    
    if (ch > '\0' ) 
      {
	std::cout << ch ; 
	--n; 	
      }
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


static boolean  treeGetLabel (std::istringstream &iss , char *lblPtr, int maxlen)
{
  int      ch;
  boolean  done, quoted, lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *) NULL; 
  else 
    if (lblPtr == NULL) 
      maxlen = 0;

  ch = iss.get();
  done = treeLabelEnd(ch);

  lblfound = ! done;
  quoted = (ch == '\'');
  if (quoted && ! done) 
    {
      ch = iss.get(); 
      done = (ch == EOF);
    }

  while (! done) 
    {
      if (quoted) 
	{
	  if (ch == '\'') 
	    {
	      ch = iss.get(); 
	      if (ch != '\'') 
		break;
	    }
        }
      else 
	if (treeLabelEnd(ch)) break;     

      if (--maxlen >= 0) *lblPtr++ = ch;
      ch = iss.get();
      if (ch == EOF) break;
    }

  if (ch != EOF)   iss.unget();

  if (lblPtr != NULL) *lblPtr = '\0';

  return lblfound;
}


static boolean  treeFlushLabel (std::istringstream &iss)
{ 
  return  treeGetLabel(iss, (char *) NULL, (int) 0);
}





static int treeFinishCom ( std::istringstream &iss, char **strp)
{
  int  ch;
  
  while ((ch = iss.get()) != EOF && ch != ']') {
    if (strp != NULL) *(*strp)++ = ch;    /* save character  */
    if (ch == '[') {                      /* nested comment; find its end */
      if ((ch = treeFinishCom(iss, strp)) == EOF)  break;
      if (strp != NULL) *(*strp)++ = ch;  /* save closing ]  */
    }
  }
  
  if (strp != NULL) **strp = '\0';        /* terminate string  */
  return  ch;
} 

static int treeGetCh ( std::istringstream &iss )         /* get next nonblank, noncomment character */
{
  int  ch;

  while ((ch = iss.get() ) != EOF) {
    if (whitechar(ch)) ;
    else if (ch == '[')
      {
	if ((ch = treeFinishCom(iss, (char **) NULL)) == EOF)  
	  break;
      }
    else  
      break;
  }
  
  return  ch;
}



static boolean treeNeedCh ( std::istringstream &iss , int c1, std::string where)
{
  int  c2;
  
  if ((c2 = treeGetCh(iss)) == c1)  return TRUE;
  
  printf("ERROR: Expecting '%c' %s tree; found:", c1, where.c_str());
  if (c2 == EOF) 
    {
      printf("End-of-File");
    }
  else 
    {      	
      iss.unget();
      treeEchoContext(iss, 40);
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


static int treeFindTipName(std::istringstream &iss , tree *tr)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabel(iss, str, nmlngth+2))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   

  return  n;
}


static boolean treeProcessLength (std::istringstream &iss, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(iss)) == EOF)  return FALSE;    /*  Skip comments */
  iss.unget();
  

  if ( iss >> *dptr ) 
    {
      return TRUE; 
    } 
  else 
    {
      printf("ERROR: treeProcessLength: Problem reading branch length\n");
      treeEchoContext(iss, 40);
      printf("\n");
      return  FALSE;
    } 

  return  TRUE;
}



static int treeFlushLen (std::istringstream &iss)
{
  double  dummy;  
  int     ch;
  
  ch = treeGetCh(iss);
  
  if (ch == ':') 
    {
      ch = treeGetCh(iss);

      iss.unget();
      if(!treeProcessLength(iss, & dummy)) return 0;
      return 1;	  
    }
  
  if (ch != EOF) 
    iss.unget();
  return 1;
}


static boolean addElementLen (std::istringstream &iss, tree *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount)
{   
  nodeptr  q;
  int      n, ch, fres;
  
  if ((ch = treeGetCh(iss)) == '(') 
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

      if (! addElementLen(iss, tr, q->next, readBranchLengths, readNodeLabels, lcount))        return FALSE;
      if (! treeNeedCh(iss, ',', "in"))             return FALSE;
      if (! addElementLen(iss, tr, q->next->next, readBranchLengths, readNodeLabels, lcount))  return FALSE;
      if (! treeNeedCh(iss, ')', "in"))             return FALSE;
      
      if(readNodeLabels)
	{
	  char label[64];
	  int support;

	  if(treeGetLabel (iss, label, 10))
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
	treeFlushLabel(iss);
    }
  else 
    {   
      iss.unget();
      if ((n = treeFindTipName(iss, tr)) <= 0)          return FALSE;
      q = tr->nodep[n];
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
    }
  
  if(readBranchLengths)
    {
      double branch;
      if (! treeNeedCh(iss, ':', "in"))                 return FALSE;
      if (! treeProcessLength(iss, &branch))            return FALSE;

      // TODO only works for 1 trbl
      hookup(p, q, &branch, 1 );
    }
  else
    {
      fres = treeFlushLen(iss);
      if(!fres) return FALSE;
      
      exa_hookupDefault(tr, p, q);
    }
  return TRUE;          
}




void myTreeReadLen(std::string treeString , tree *tr, boolean hasBL)
{
  auto&& iss = std::istringstream{treeString}; 

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
  // for(i = 0; i < 1; i++)
  //   tr->partitionSmoothed[i] = FALSE;
  
  tr->rooted      = FALSE;     

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(iss)) != '(');


  /* if(!topologyOnly) */
  /*   assert(readBranches == FALSE && readNodeLabels == FALSE); */
  
       
  if (! addElementLen(iss, tr, p, readBranches, readNodeLabels, &lcount))                 
    assert(0);
  if (! treeNeedCh(iss, ',', "in"))                
    assert(0);
  if (! addElementLen(iss, tr, p->next, readBranches, readNodeLabels, &lcount))
    assert(0);
  if (! tr->rooted) 
    {
      if ((ch = treeGetCh(iss)) == ',') 
	{ 
	  if (! addElementLen(iss, tr, p->next->next, readBranches, readNodeLabels, &lcount))
	    assert(0);	    
	}
      else 
	{                                    /*  A rooted format */
	  tr->rooted = TRUE;
	  if (ch != EOF) 
	    iss.unget();
	}	
    }
  else 
    {      
      p->next->next->back = (nodeptr) NULL;
    }
  if (! treeNeedCh(iss, ')', "in"))                
    assert(0);

  if(topologyOnly)
    assert(!(tr->rooted && readNodeLabels));

  treeFlushLabel(iss);
  
  if (! treeFlushLen(iss))                         
    assert(0);
 
  if (! treeNeedCh(iss, ';', "at end of"))       
    assert(0);

  tr->start = findAnyTip(p, tr->mxtips);    
 
  assert(tr->ntips == tr->mxtips);
}


boolean readTreeWithOrWithoutBL(tree *tr, std::string treeString)
{  
  bool hasBL = std::any_of(treeString.begin(), treeString.end(), [](char c){ return c == ':' ; }); 
  myTreeReadLen(treeString, tr, hasBL); 

  return hasBL; 
}
