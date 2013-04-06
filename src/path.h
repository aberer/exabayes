#ifndef _PATH_H
#define _PATH_H


#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "adapters.h"
#include "branch.h"


typedef stack path;  

void pushToStackIfNovel(stack *s, branch b, int numTip); 
void drawPathForESPR(state *chain,stack *s, double stopProp ); 
void saveBranchLengthsPath(state *chain, path *s); 
boolean nodeIsOnPath(int node, path *aPath); 
boolean isOuterNode(int node, path *aPath); 
void multiplyAlongBranchESPR(state *chain, path *s, double multi); 
#endif
