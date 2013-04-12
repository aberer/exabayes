/**
   @file exa-topology.h
    
   @brief Adapted code from the axml-variants for storing topologies.  


   TODO on the long run this should be replaced by a version tailored
   to our needs. The axml-topology does not work across trees and
   therefore many hacks in this file could be omitted.
   
*/

#ifndef _BAYES_TOPO
#define _BAYES_TOPO


/* #ifdef __cplusplus */
/* extern "C"{ */
/* #endif */


topol  *setupTopol (int maxtips);
void  freeTopol (topol *tpl);
void saveTree (tree *tr, topol *tpl); 
boolean restoreTree (topol *tpl, tree *tr);
void copyTopology(tree *targetTree, tree *origTree); 


/* #ifdef __cplusplus */
/* } */
/* #endif */

#endif
