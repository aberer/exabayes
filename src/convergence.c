#include "globals.h"
#include "common.h"

#include "config.h"
#include "axml.h"

#include "proposalStructs.h"
#include "main-common.h"






/* TODO does not diagnose anything currently    */
boolean convergenceDiagnosticDummy(state *allChains, int numChains)
{
  boolean  result = TRUE; 
  
  for(int i = 0; i < numChains; ++i)
    {
      state *curChain = allChains + i; 
      if(curChain->currentGeneration < curChain->numGen)
	result = FALSE; 
    }
  
  return result; 
}



boolean convergenceDiagnostic(state *allChains, int numChains)
{
  return convergenceDiagnosticDummy(allChains, numChains); 
}
