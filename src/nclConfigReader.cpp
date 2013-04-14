

#include <iostream>
#include <cassert>
#include <ncl/ncl.h>

#include "nclConfigReader.h"


// TODO set this up more generic   


class ConfigReader : public NxsReader
{
public: 
  ConfigReader() : NxsReader(){SetWarningOutputLevel(SUPPRESS_WARNINGS_LEVEL); }
  virtual void ExitingBlock(NxsString blockName){}
  virtual void ExecuteStopping(){}
  virtual void ExecuteStarting(){}
  
} ; 


/**
   @brief This function is nothing to be proud of, but let's stick to
   it for now.
 */ 
bool mapNameToProposal(NxsString &key, proposal_type *pf)
{    
  bool found = true; 

  if(key.EqualsCaseInsensitive("initModelWeight")) 
    *pf = UPDATE_MODEL; 
  else if(key.EqualsCaseInsensitive("initGammaWeight")) 
    *pf = UPDATE_GAMMA; 
  else if(key.EqualsCaseInsensitive("treeLengthMult")) 
    *pf = TL_MULT; 
  else if(key.EqualsCaseInsensitive("rateHetMulti"))    
    *pf = GAMMA_MULTI; 
  else if(key.EqualsCaseInsensitive("initGammaExpWeight"))    
    *pf = UPDATE_GAMMA_EXP; 
  else if(key.EqualsCaseInsensitive("initSingleBranchWeight"))    
    *pf = UPDATE_SINGLE_BL; 
  else if(key.EqualsCaseInsensitive("initSingleBranchExpWeight"))    
    *pf = UPDATE_SINGLE_BL_EXP; 
  else if(key.EqualsCaseInsensitive("initSingleBranchBiunifWeight"))    
    *pf = UPDATE_SINGLE_BL_BIUNIF; 
  else if(key.EqualsCaseInsensitive("initModelBiunifWeight"))    
    *pf = UPDATE_MODEL_BIUNIF; 
  else if(key.EqualsCaseInsensitive("initModelSingleBiunifWeight"))    
    *pf = UPDATE_MODEL_SINGLE_BIUNIF; 
  else if(key.EqualsCaseInsensitive("initModelAllBiunifWeight"))    
    *pf = UPDATE_MODEL_ALL_BIUNIF; 
  else if(key.EqualsCaseInsensitive("initModelPermBiunifWeight"))    
    *pf = UPDATE_MODEL_PERM_BIUNIF; 
  else if(key.EqualsCaseInsensitive("initFrequenciesWeight"))    
    *pf = UPDATE_FREQUENCIES_BIUNIF; 
  else if(key.EqualsCaseInsensitive("initEsprMappedWeight"    ))    
    *pf = E_SPR; 
  else if(key.EqualsCaseInsensitive("branchMulti"))
    *pf = BRANCH_LENGTHS_MULTIPLIER; 
  else if(key.EqualsCaseInsensitive("initFrequencySliderWeight"))
    *pf = FREQUENCY_SLIDER; 
  else if(key.EqualsCaseInsensitive("initGuidedSPR"))
    *pf = GUIDED_SPR; 
  else if(key.EqualsCaseInsensitive("stNNI"))
    *pf = ST_NNI; 
  else if (key.EqualsCaseInsensitive("nodeSlider"))
    *pf = NODE_SLIDER; 
  // TODO@kassian this is a good place for proposal add 

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
    initParam->initPenaltyFactor  = -1 ;   
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
	    //PROPOSALADD read NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
	    else if(key.EqualsCaseInsensitive("numGen"))
	      initParam->numGen = value.ConvertToInt(); 
	    else if(key.EqualsCaseInsensitive("initPenaltyFactor"))
	      initParam->initPenaltyFactor = value.ConvertToDouble();
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
    assert(initParam->initPenaltyFactor != -1); 
    assert(initParam->samplingFrequency  != -1 ); 
    assert(initParam->eSprStopProb != -1); 
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



void parseConfigWithNcl(char *configFileName, initParamStruct **params)
{
  ConfigReader reader;
  ExabayesBlock *myBlock = new ExabayesBlock(); 

  ifstream fh(configFileName); 
  
  NxsToken token(fh); 

  reader.Add(myBlock); 
  reader.Execute(token);
  
  *params = myBlock->getResult();
}

