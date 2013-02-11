#include <iostream>
#include <cassert>

#include <ncl/ncl.h>
#include "nclConfigReader.h"

// NOTE: this is a c++ file! 


class ConfigReader : public NxsReader
{
public: 
  ConfigReader() : NxsReader(){SetWarningOutputLevel(SUPPRESS_WARNINGS_LEVEL); }
  virtual void ExitingBlock(NxsString blockName){}
  virtual void ExecuteStopping(){}
  virtual void ExecuteStarting(){}
  
} ; 



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
    initParam->initSPRWeight  = -1 ;
    initParam->initGammaWeight  = -1 ; 
    initParam->initModelWeight  = -1 ;
    initParam->initSingleBranchWeight  = -1 ; 
    initParam->initSingleBranchExpWeight  = -1 ;
    initParam->initPenaltyFactor  = -1 ;   
    initParam->numGen  = -1 ; 
    initParam->samplingFrequency  = -1 ; 
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

	    if( key.EqualsCaseInsensitive("initSPRWeight") )
	      initParam->initSPRWeight = value.ConvertToDouble();
	    else if(key.EqualsCaseInsensitive("initModelWeight"))
	      initParam->initModelWeight = value.ConvertToDouble(); 
	    else if(key.EqualsCaseInsensitive("initGammaWeight"))
	      initParam->initGammaWeight = value.ConvertToDouble();
	    else if (key.EqualsCaseInsensitive("initSingleBranchWeight"))
	      initParam->initSingleBranchWeight = value.ConvertToDouble();
	    else if (key.EqualsCaseInsensitive("initSingleBranchExpWeight"))
	      initParam->initSingleBranchExpWeight = value.ConvertToDouble();	    
	    else if(key.EqualsCaseInsensitive("numGen"))
	      initParam->numGen = value.ConvertToInt(); 
	    else if(key.EqualsCaseInsensitive("initPenaltyFactor"))
	      initParam->initPenaltyFactor = value.ConvertToDouble();
	    else if(key.EqualsCaseInsensitive("samplingfrequency"))
	      initParam->samplingFrequency = value.ConvertToInt();
	    else 	      
	      cerr << "WARNING: ignoring unknown value >"  << key << "< and >" << value <<  "<" << endl; 
	  }
      }
  }

  void assertInitialized()
  {
    assert(initParam->initModelWeight != -1); 
    assert(initParam->initGammaWeight != -1); 
    assert(initParam->initSPRWeight != -1); 
    assert(initParam->initSingleBranchWeight != -1); 
    assert(initParam->numGen != -1); 
    assert(initParam->initPenaltyFactor != -1); 
    assert(initParam->samplingFrequency  != -1 ); 
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

