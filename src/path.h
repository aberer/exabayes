#ifndef _PATH_H
#define _PATH_H


#include "axml.h"
#include "bayes.h"
/* #include "randomness.h" */
#include "adapters.h"
#include "branch.h"


typedef stack path;  

void pushToStackIfNovel(stack *s, branch b, int numTip); 
void drawPathForESPR(Chain *chain,stack *s, double stopProp ); 
void saveBranchLengthsPath(Chain *chain, path *s); 
boolean nodeIsOnPath(int node, path *aPath); 
boolean isOuterNode(int node, path *aPath); 
void multiplyAlongBranchESPR(Chain *chain, path *s, double multi); 
#endif
