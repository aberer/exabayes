
#pragma   once 


#include "axml.h"
#include "bayes.h"

#include "randomness.h"
#include "TreeAln.hpp"

void applyNodeSlider(state *chain, proposalFunction *pf); 
void evaluateNodeSlider(state *chain, proposalFunction *pf); 
void resetNodeSlider(state *chain, proposalFunction *pf); 
