#include "BlockProposalConfig.hpp"

#include <cassert>


BlockProposalConfig::BlockProposalConfig()
  : setByUser(NUM_PROPOSALS, false)
   ,userProposalWeights(NUM_PROPOSALS, 0)
{
  NCL_BLOCKTYPE_ATTR_NAME = "PROPOSALS"; 
  setupMap();
}


void BlockProposalConfig::Read(NxsToken &token)
{ 
  DemandEndSemicolon(token, "PROPOSALS");

  while(true)
    {
      token.GetNextToken();
      NxsBlock::NxsCommandResult res = HandleBasicBlockCommands(token); 
   
      if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
	return;
      if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND))
	{
	  NxsString key = token.GetToken(false);
	  token.GetNextToken(); 
	  NxsString value = token.GetToken(false); 	    

	  if(name2proposal.find(key) != name2proposal.end() )
	    {
	      double val = value.ConvertToDouble(); 
	      userProposalWeights[proposal_type(name2proposal[key])] = val;
	      setByUser[proposal_type(name2proposal[key])] = true; 
	    }
	  else if(key.EqualsCaseInsensitive("esprstopprob"))	    
	    esprStopProp = value.ConvertToDouble();	  
	  else if(key.EqualsCaseInsensitive("guidedsprradius"))
	    guidedRadius = value.ConvertToInt();
	  else if(key.EqualsCaseInsensitive("parsimonyWarp"))	    
	    parsimonyWarp = value.ConvertToDouble();
	  else 	      
	    {
	      cerr << "WARNING: ignoring unknown value >"  << key << "< and >" << value <<  "<" << endl; 
	      assert(0);
	    }
	}
    }
}


// NOTICE 
void BlockProposalConfig::setupMap()
{

  // TODO 
  name2proposal["STNNI"] =  ST_NNI;
  name2proposal["ESPR"] =  E_SPR;
  name2proposal["ETBR"] =  E_TBR;
  name2proposal["PARSIMONYSPR"] =  PARSIMONY_SPR;
  name2proposal["GUIDEDSPR"] =  GUIDED_SPR;
  
  // BL 
  name2proposal["BRANCHSLIDER"] =  BRANCH_SLIDER;
  name2proposal["BRANCHCOLLAPSER"] = BRANCH_COLLAPSER; 
  name2proposal["TREELENGTHMULT"] =  TL_MULT;
  name2proposal["BRANCHMULTI"] =  BRANCH_LENGTHS_MULTIPLIER;
  name2proposal["GUIDEDBL"] =  UPDATE_SINGLE_BL_GUIDED;
  name2proposal["NODESLIDER"] =  NODE_SLIDER;
  
  // revmat 
  name2proposal["REVMATSLIDER"] = REVMAT_SLIDER; 
  name2proposal["REVMATDIRICHLET"] =  REVMAT_DIRICHLET;
  
  // rate heterogeneity 
  name2proposal["RATEHETSLIDER"] =  RATE_HET_SLIDER;
  name2proposal["RATEHETMULTI"] =  RATE_HET_MULTI;

  // state frequencies
  name2proposal["FREQUENCYSLIDER"] =  FREQUENCY_SLIDER;
  name2proposal["FREQUENCYDIRICHLET"] =  FREQUENCY_DIRICHLET;

  name2proposal["AAMODELJUMP"] = AMINO_MODEL_JUMP; 
}
