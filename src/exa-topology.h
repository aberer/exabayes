/**
   @file exa-topology.h
    
   @brief Adapted code from the axml-variants for storing topologies.  


   TODO on the long run this should be replaced by a version tailored
   to our needs. The axml-topology does not work across trees and
   therefore many hacks in this file could be omitted.
   
*/

#ifndef _BAYES_TOPO
#define _BAYES_TOPO


topol  *setupTopol (int maxtips);
void  freeTopol (topol *tpl);
void saveTree (TreeAln *traln, topol *tpl); 
boolean restoreTree (topol *tpl, TreeAln *traln);
void copyTopology(tree *targetTree, tree *origTree); 

#endif
