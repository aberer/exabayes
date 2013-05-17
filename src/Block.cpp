#include "Block.hpp"


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
, proposalWeights(NUM_PROPOSALS)
{
  NCL_BLOCKTYPE_ATTR_NAME = "EXABAYES"; 
  setupMap();
}


shared_ptr<AbstractPrior> ExabayesBlock::parsePrior(NxsToken &token, NxsString &value)
{
  token.GetNextToken();
  assert(token.GetToken().compare("(") == 0); 
  if(value.EqualsCaseInsensitive("uniform"))
    {
      token.GetNextToken();

      // for non-continuous variables (e.g., topology)
      if(token.GetToken().compare(")") == 0) 
	return shared_ptr<AbstractPrior> (new UniformPrior(0,0)); 

      double n1 = atof(token.GetToken().c_str()); 
      token.GetNextToken();
      assert(token.GetToken().compare(",") == 0);
      token.GetNextToken();
      double n2 = atof(token.GetToken().c_str());
      token.GetNextToken();
      assert(token.GetToken().compare(")") == 0);
      return shared_ptr<AbstractPrior> (new UniformPrior(n1,n2));  
    }
  else if(value.EqualsCaseInsensitive("dirichlet"))
    {      
      vector<double> alphas; 
      while(token.GetToken().compare(")") != 0 )
	{
	  token.GetNextToken();
	  alphas.push_back(atof(token.GetToken().c_str())); 
	  token.GetNextToken();
	  assert(token.GetToken().compare(",") == 0 
		 || token.GetToken().compare(")") == 0); 
	} 
      return shared_ptr<AbstractPrior>(new DirichletPrior(alphas)); 
    }
  else if(value.EqualsCaseInsensitive("fixed"))
    { 
      token.GetNextToken();
      if(token.GetToken().EqualsCaseInsensitive("empirical"))
	{
	  cerr << "not implemented yet " << endl; // TODO 
	  exit(1);
	  return NULL; 
	} 
      else 
	{
	  vector<double> fixedValues; 	  
	  while(token.GetToken().compare(")") != 0)
	    {
	      fixedValues.push_back(atof(token.GetToken().c_str()));
	      token.GetNextToken();
	      if(token.GetToken().compare(",") == 0)
		token.GetNextToken();
	    }
	  return shared_ptr<AbstractPrior>(new FixedPrior(fixedValues));
	}
    }
  else if(value.EqualsCaseInsensitive("exponential"))
    {
      token.GetNextToken();
      double n1 = atof(token.GetToken().c_str());
      token.GetNextToken();
      assert(token.GetToken().compare(")") == 0);
      return shared_ptr<AbstractPrior>(new ExponentialPrior(n1));
    }
  else 
    {
      cerr << "attempted to parse prior. Did not recognize keyword " <<  value << endl; 
      exit(1);
    }

  return NULL; 
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
	  else if(key.EqualsCaseInsensitive("esprstopprob"))	    
	    sm->setEsprStopProp( value.ConvertToDouble());
	  else if(key.EqualsCaseInsensitive("numIndiChains"))
	    sm->setNumRunConv(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("diagFreq"))
	    sm->setDiagFreq(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("numcoupledChains"))
	    sm->setNumCoupledChains(value.ConvertToInt());
	  else if(key.EqualsCaseInsensitive("guidedsprradius"))
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
	  // PRIORS
	  else if(key.EqualsCaseInsensitive("brlenpr"))
	    prior.setBranchLengthPrior(parsePrior(token, value)); 
	  else if(key.EqualsCaseInsensitive("revmatpr"))
	    prior.setRevMatPrior( parsePrior(token,value)); 
	  else if(key.EqualsCaseInsensitive("statefreqpr"))
	    prior.setStateFreqPrior(parsePrior(token,value ));
	  else if(key.EqualsCaseInsensitive("shapepr"))	      
	    prior.setRateHetPrior(parsePrior(token,value)); 	    
	  else if(key.EqualsCaseInsensitive("topoPr"))
	    prior.setTopologyPrior(parsePrior(token,value)); 
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
}

