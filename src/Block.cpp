#include "Block.hpp"

#include "proposalType.h"


static bool convertToBool(NxsString &string)
{
  if(string.EqualsCaseInsensitive("true"))
    return true ; 
  else if (string.EqualsCaseInsensitive("false"))
    return false; 
  else 
    {
      cerr << "ERROR while parsing boolean value: expected either \"true\" or \"false\"" << endl; 
      assert(0); 
      return false; 
    }  
}


ExabayesBlock::ExabayesBlock(SampleMaster *_sm ) 
:sm(_sm)
, proposalWeights(0, NUM_PROPOSALS)
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
	  else if(key.EqualsCaseInsensitive("numGen"))
	    sm->setNumGen(value.ConvertToInt()); 
	  else if(key.EqualsCaseInsensitive("samplingfrequency"))
	    sm->setSamplingFreq(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("eSprStopProb"))
	    sm->setEsprStopProp( value.ConvertToDouble());
	  else if(key.EqualsCaseInsensitive("numIndiChains"))
	    sm->setNumRunConv(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("diagFreq"))
	    sm->setDiagFreq(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("numcoupledChains"))
	    sm->setNumCoupledChains(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("guidedSPRRadius"))
	    sm->setGuidedRadius(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("printFreq") )	   
	    sm->setPrintFreq( value.ConvertToInt());
	  else if (key.EqualsCaseInsensitive("asdsfIgnoreFreq"))
	    sm->setAsdsfIgoreFreq( value.ConvertToDouble());
	  else if (key.EqualsCaseInsensitive("asdsfConvergence"))
	    sm->setAsdsfConvergence( value.ConvertToDouble());
	  else if (key.EqualsCaseInsensitive("heatFactor"))
	    sm->setHeatFactor(value.ConvertToDouble());
	  else if(key.EqualsCaseInsensitive("swapInterval"))
	    sm->setSwapInterval( value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("tuneHeat"))
	    sm->setTuneHeat( convertToBool(value));
	  else if(key.EqualsCaseInsensitive("tuneFreq"))
	    sm->setTuneFreq(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("burninGen"))
	    sm->setBurninGen( value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("burninProportion"))
	    sm->setBurninProportion( value.ConvertToDouble());
	  else if(key.EqualsCaseInsensitive("parsimonyWarp"))
	    sm->setParsimonyWarp( value.ConvertToDouble());
	  else if(key.EqualsCaseInsensitive("brlen"))
	    {
	      cout << "yes! value would be "<< value << endl; 

	      token.GetNextToken();		
	      string tok = token.GetToken();
	      assert(tok.compare ("(") == 0 );
	      token.GetNextToken();		
	      // tok = 
	    }
	  else 	      
	    cerr << "WARNING: ignoring unknown value >"  << key << "< and >" << value <<  "<" << endl; 
	}
    }
}





void ExabayesBlock::setupMap()
{
  name2proposal["model"] = UPDATE_MODEL; 
  name2proposal["ratehetslider"] =  UPDATE_GAMMA;
  name2proposal["treelengthmult"] =  TL_MULT;
  name2proposal["ratehetmulti"] =  GAMMA_MULTI;
  name2proposal["ratehetexp"] =  UPDATE_GAMMA_EXP;
  name2proposal["singlebranch"] =  UPDATE_SINGLE_BL;
  name2proposal["singlebranchexp"] =  UPDATE_SINGLE_BL_EXP;
  name2proposal["singlebranchbiunif"] =  UPDATE_SINGLE_BL_BIUNIF;
  name2proposal["modelbiunif"] =  UPDATE_MODEL_BIUNIF;
  name2proposal["modelsinglebiunif"] =  UPDATE_MODEL_SINGLE_BIUNIF;
  name2proposal["modelallbiunif"] =  UPDATE_MODEL_ALL_BIUNIF;
  name2proposal["modelpermbiunif"] =  UPDATE_MODEL_PERM_BIUNIF;
  name2proposal["modeldirichlet"] =  UPDATE_MODEL_DIRICHLET;
  name2proposal["frequencies"] =  UPDATE_FREQUENCIES_BIUNIF;
  name2proposal["espr"] =  E_SPR;
  name2proposal["branchmulti"] =  BRANCH_LENGTHS_MULTIPLIER;
  name2proposal["guidedbl"] =  UPDATE_SINGLE_BL_GUIDED;
  name2proposal["frequencyslider"] =  FREQUENCY_SLIDER;
  name2proposal["guidedspr"] =  GUIDED_SPR;
  name2proposal["stnni"] =  ST_NNI;
  name2proposal["nodeslider"] =  NODE_SLIDER;
  name2proposal["frequencydirichlet"] =  UPDATE_FREQUENCIES_DIRICHLET;
  name2proposal["etbr"] =  E_TBR;
  name2proposal["parsimonyspr"] =  PARSIMONY_SPR;
}

