#ifndef _BAYES_TOPO
#define _BAYES_TOPO


topol  *setupTopol (int maxtips);
void  freeTopol (topol *tpl);
void saveTree (tree *tr, topol *tpl); 
boolean restoreTree (topol *tpl, tree *tr);
/* void copyTopoFromDifferentTree(tree *targetTree, tree *otherTree, nodeptr p) ;  */
void copyTopology(tree *targetTree, tree *origTree); 

#endif
