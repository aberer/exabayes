#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "main-common.h"
#include "TreeAln.hpp" 
#include "adapters.h"
#include "BipartitionHash.hpp"

void printBv(unsigned int *bv, int bvlen)
{
  if(isOutputProcess())
    {
      for(int i = 0; i < bvlen; ++i)
	{
	  int pos = i / 32; 
	  int relPos = i % 32; 
	if( ( bv[pos ] & (1 << (relPos)) ) > 0 )
	  printf("1"); 
	else 
	  printf("0");
	}
    }
}


/* TODO does not diagnose anything currently    */
boolean convergenceDiagnosticDummy(state *allChains, int numChains)
{
  boolean  result = TRUE; 
  
  for(int i = 0; i < numChains; ++i)
    {
      state *curChain = allChains + i; 
      if(curChain->currentGeneration < gAInfo.numGen)
	result = FALSE; 
    }
  
  return result; 
}


boolean convergenceDiagnostic(state *allChains)
{
  if(gAInfo.numberOfRuns > 1)
    {
      double asdsf = gAInfo.bipHash->averageDeviationOfSplitFrequencies(gAInfo.asdsfIgnoreFreq);
      PRINT("ASDSF=%f\n", asdsf); 
      return asdsf < gAInfo.asdsfConvergence ; 
    }
  else 
    return allChains[0].currentGeneration > gAInfo.numGen; 
}



