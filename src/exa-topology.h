#ifndef _BAYES_TOPO
#define _BAYES_TOPO


void traverseAndCount(nodeptr p, int *count, tree *tr ); 
topol  *setupTopol (int maxtips);
void  freeTopol (topol *tpl);
void saveTree (tree *tr, topol *tpl); 
boolean restoreTree (topol *tpl, tree *tr);
void copyTopology(tree *targetTree, tree *origTree); 

#endif
