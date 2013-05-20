#include "Block.hpp"


ExabayesBlock::ExabayesBlock(SampleMaster *_sm ) 
:sm(_sm)
, proposalWeights(NUM_PROPOSALS)
{
  NCL_BLOCKTYPE_ATTR_NAME = "EXABAYES"; 
  setupMap();
}






void ExabayesBlock::Read(NxsToken &token)
{ 
  DemandEndSemicolon(token, "EXABAYES");

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
	    proposalWeights[proposal_type(name2proposal[key])] = value.ConvertToDouble(); 
	  // else if(key.EqualsCaseInsensitive("numGen"))
	    // sm->setNumGen(value.ConvertToInt()); 
	  // else if(key.EqualsCaseInsensitive("samplingfrequency"))
	    // sm->setSamplingFreq(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("esprstopprob"))	    
	    sm->setEsprStopProp( value.ConvertToDouble());
	  // else if(key.EqualsCaseInsensitive("numIndiChains"))
	  //   sm->setNumRunConv(value.ConvertToInt());
	  // else if(key.EqualsCaseInsensitive("diagFreq"))
	  //   sm->setDiagFreq(value.ConvertToInt());
	  // else if(key.EqualsCaseInsensitive("numcoupledChains"))
	  //   sm->setNumCoupledChains(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("guidedsprradius"))
	    sm->setGuidedRadius(value.ConvertToInt());
	  // else if(key.EqualsCaseInsensitive("printFreq") )	   
	  //   sm->setPrintFreq( value.ConvertToInt());
	  // else if (key.EqualsCaseInsensitive("asdsfIgnoreFreq"))
	  //   sm->setAsdsfIgoreFreq( value.ConvertToDouble());
	  // else if (key.EqualsCaseInsensitive("asdsfConvergence"))
	  //   sm->setAsdsfConvergence( value.ConvertToDouble());
	  // else if (key.EqualsCaseInsensitive("heatFactor"))
	  //   sm->setHeatFactor(value.ConvertToDouble());
	  // else if(key.EqualsCaseInsensitive("swapInterval"))
	  //   sm->setSwapInterval( value.ConvertToInt());
	  // else if(key.EqualsCaseInsensitive("tuneHeat"))
	  //   sm->setTuneHeat( convertToBool(value));
	  // else if(key.EqualsCaseInsensitive("tuneFreq"))
	  //   sm->setTuneFreq(value.ConvertToInt());
	  // else if(key.EqualsCaseInsensitive("burninGen"))
	  //   sm->setBurninGen( value.ConvertToInt());
	  // else if(key.EqualsCaseInsensitive("burninProportion"))
	  //   sm->setBurninProportion( value.ConvertToDouble());
	  else if(key.EqualsCaseInsensitive("parsimonyWarp"))
	    sm->setParsimonyWarp( value.ConvertToDouble());
	  else 	      
	    cerr << "WARNING: ignoring unknown value >"  << key << "< and >" << value <<  "<" << endl; 
	}
    }
}




// NOTICE 
void ExabayesBlock::setupMap()
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

