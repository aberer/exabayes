
#include "config.h"

#include <iostream>
#include <cassert>
#include "ConfigReader.hpp"


// TODO set this up more generic   


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


/**
   @brief This function is nothing to be proud of, but let's stick to
   it for now.
 */ 
bool mapNameToProposal(NxsString &key, proposal_type *pf)
{    
  bool found = true; 

  if(key.EqualsCaseInsensitive("model")) 
    *pf = UPDATE_MODEL; 
  else if(key.EqualsCaseInsensitive("rateHetslider")) 
    *pf = UPDATE_GAMMA; 
  else if(key.EqualsCaseInsensitive("treeLengthMult")) 
    *pf = TL_MULT; 
  else if(key.EqualsCaseInsensitive("rateHetMulti"))    
    *pf = GAMMA_MULTI; 
  else if(key.EqualsCaseInsensitive("ratehetexp"))    
    *pf = UPDATE_GAMMA_EXP; 
  else if(key.EqualsCaseInsensitive("singlebranch"))    
    *pf = UPDATE_SINGLE_BL; 
  else if(key.EqualsCaseInsensitive("singlebranchexp"))    
    *pf = UPDATE_SINGLE_BL_EXP; 
  else if(key.EqualsCaseInsensitive("singlebranchbiunif"))    
    *pf = UPDATE_SINGLE_BL_BIUNIF; 
  else if(key.EqualsCaseInsensitive("modelbiunif"))    
    *pf = UPDATE_MODEL_BIUNIF; 
  else if(key.EqualsCaseInsensitive("modelsinglebiunif"))    
    *pf = UPDATE_MODEL_SINGLE_BIUNIF; 
  else if(key.EqualsCaseInsensitive("modelallbiunif"))    
    *pf = UPDATE_MODEL_ALL_BIUNIF; 
  else if(key.EqualsCaseInsensitive("modelpermbiunif"))    
    *pf = UPDATE_MODEL_PERM_BIUNIF;
  else if(key.EqualsCaseInsensitive("modeldirichlet"))    
    *pf = UPDATE_MODEL_DIRICHLET; 
  else if(key.EqualsCaseInsensitive("frequencies"))    
    *pf = UPDATE_FREQUENCIES_BIUNIF; 
  else if(key.EqualsCaseInsensitive("eSPR" ))    
    *pf = E_SPR; 
  else if(key.EqualsCaseInsensitive("branchMulti"))
    *pf = BRANCH_LENGTHS_MULTIPLIER; 
  else if(key.EqualsCaseInsensitive("guidedbl"))
    *pf = UPDATE_SINGLE_BL_GUIDED; 
  else if(key.EqualsCaseInsensitive("frequencyslider"))
    *pf = FREQUENCY_SLIDER; 
  else if(key.EqualsCaseInsensitive("guidedSPR"))
    *pf = GUIDED_SPR; 
  else if(key.EqualsCaseInsensitive("stNNI"))
    *pf = ST_NNI; 
  else if (key.EqualsCaseInsensitive("nodeSlider"))
    *pf = NODE_SLIDER; 
  else if (key.EqualsCaseInsensitive("frequencydirichlet"))
    *pf = UPDATE_FREQUENCIES_DIRICHLET;
  else if (key.EqualsCaseInsensitive("etbr"))
    *pf = E_TBR;
  else if(key.EqualsCaseInsensitive("parsimonySPR"))
    *pf = PARSIMONY_SPR; 
  else 
    found = false; 

  return found; 
}


// TODO: catch exceptions: for now, we assume that our users are
// capable of specifying the file correctly
class ExabayesBlock : public NxsBlock
{
public: 
  ExabayesBlock(void ) 
  {
    NCL_BLOCKTYPE_ATTR_NAME = "EXABAYES"; 
    initParam = (initParamStruct*)calloc(1,sizeof(initParamStruct)); 

    // note: all values are -1 by default, so we can check for meaningful values later 
    paramBadInit();
  }

  void paramBadInit()
  {

    initParam->numCoupledChains = -1; 
    initParam->numIndiChains = -1; 
    initParam->numGen  = -1 ; 
    initParam->samplingFrequency  = -1 ; 
    initParam->eSprStopProb = -1;
    initParam->diagFreq = -1; 

  }


  virtual void Read(NxsToken &token)
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
	    proposal_type pt; 
	    
	    if(mapNameToProposal(key, &pt))
	      initParam->initWeights[pt] = value.ConvertToDouble(); 
	    else if(key.EqualsCaseInsensitive("numGen"))
	      initParam->numGen = value.ConvertToInt(); 
	    else if(key.EqualsCaseInsensitive("samplingfrequency"))
	      initParam->samplingFrequency = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("eSprStopProb"))
	      initParam->eSprStopProb = value.ConvertToDouble();
	    else if(key.EqualsCaseInsensitive("numIndiChains"))
	      initParam->numIndiChains = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("diagFreq"))
	      initParam->diagFreq = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("numCoupledChains"))
	      initParam->numCoupledChains = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("guidedSPRRadius"))
	      initParam->initGuidedSPR = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("printFreq"))
	      initParam->printFreq = value.ConvertToInt();
	    else if (key.EqualsCaseInsensitive("asdsfIgnoreFreq"))
	      initParam->asdsfIgnoreFreq = value.ConvertToDouble();
	    else if (key.EqualsCaseInsensitive("asdsfConvergence"))
	      initParam->asdsfConvergence = value.ConvertToDouble();
	    else if (key.EqualsCaseInsensitive("heatFactor"))
	      initParam->heatFactor = value.ConvertToDouble();
	    else if(key.EqualsCaseInsensitive("swapInterval"))
	      initParam->swapInterval = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("tuneHeat"))
	      initParam->tuneHeat = convertToBool(value);
	    else if(key.EqualsCaseInsensitive("tuneFreq"))
	      initParam->tuneFreq = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("burninGen"))
	      initParam->burninGen = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("burninProportion"))
	      initParam->burninProportion = value.ConvertToDouble();
	    else if(key.EqualsCaseInsensitive("numRunParallel"))
	      initParam->numRunParallel = value.ConvertToInt();
	    else if(key.EqualsCaseInsensitive("parsimonyWarp"))
	      initParam->parsWarp = value.ConvertToDouble();
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

  void assertInitialized()
  {
    assert(initParam->numCoupledChains != -1); 
    assert(initParam->diagFreq != -1 ); 
    assert(initParam->numIndiChains != -1); 
    assert(initParam->numGen != -1); 
    assert(initParam->samplingFrequency  != -1 ); 
    assert(initParam->eSprStopProb != -1); 
  }


  void validate()
  {
    if(initParam->burninProportion != 0 && initParam->burninGen)
      {
	cout << "Please EITHER specifiy the a relative burnin via burninGen OR or a relative burnin that discards a certain percentage of samples via burninProportion. Both options cannot be used at the same time." << endl; 
	exit(1); 
      }
    
    if( initParam->burninProportion > 0.99 )
      {
	cout << "Relative burnin proportion " << initParam->burninProportion << " is too high! Please choose a value below 0.99" <<endl; 
	exit(1); 
      }

    if(initParam->numIndiChains < initParam->numRunParallel)
      {
	cout << "According to your specification, " << PROGRAM_NAME << " should run " <<  initParam->numRunParallel << " runs in parallel. However, only "   << initParam->numIndiChains << " runs were specified at all. Please correct." << endl; 
	exit(1); 
      }    

  }

  initParamStruct* getResult() 
  {
    assertInitialized();
    return initParam; 
  }


private: 
  initParamStruct *initParam; 
}; 

using namespace std; 



// void parseConfigWithNcl(char *configFileName, initParamStruct **params)
// {
//   ConfigReader reader;
//   ExabayesBlock *myBlock = new ExabayesBlock(); 

//   ifstream fh(configFileName); 
  
//   NxsToken token(fh); 

//   reader.Add(myBlock); 
//   reader.Execute(token);

//   myBlock->validate();

//   *params = myBlock->getResult();
// }

