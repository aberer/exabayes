#include "BlockRunParameters.hpp" 


BlockRunParameters::BlockRunParameters()  
  : diagFreq(1000) 
  , asdsfIgnoreFreq(0.1)
  , asdsfConvergence (0.001)
  , burninGen(0)
  , burninProportion(0.25)
  , samplingFreq (100)
  , numRunConv(1)
  , numGen(50000)
  // , runId("standardId")
  , numCoupledChains(1)
  , printFreq (500)
  , heatFactor(0.1)
  , swapInterval (1)
  , tuneHeat (false)
  , tuneFreq (100)
  , useParsimonyStarting(false)
  , heatedChainsUseSame(false)
  , chkpntFreq(1000)
{
  NCL_BLOCKTYPE_ATTR_NAME = "runconfig"; 
}


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


void BlockRunParameters::Read(NxsToken &token)
{ 
  DemandEndSemicolon(token, "runconfig");

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

	  if(key.EqualsCaseInsensitive("numGen"))
	    numGen = value.ConvertToInt(); 
	  else if (key.EqualsCaseInsensitive("parsimonyStartingTree"))
	    useParsimonyStarting = convertToBool(value); 
	  else if (key.EqualsCaseInsensitive("checkpointinterval"))
	    chkpntFreq = value.ConvertToInt(); 
	  else if(key.EqualsCaseInsensitive("samplingfrequency"))
	    samplingFreq = value.ConvertToInt(); 
	  else if(key.EqualsCaseInsensitive("numRuns"))
	    numRunConv = value.ConvertToInt(); 
	  else if(key.EqualsCaseInsensitive("diagFreq"))
	    diagFreq = value.ConvertToInt(); 
	  else if(key.EqualsCaseInsensitive("heatedChainsUseSame"))
	    heatedChainsUseSame = convertToBool(value); 
	  else if(key.EqualsCaseInsensitive("numcoupledChains"))
	    numCoupledChains = value.ConvertToInt(); 
	  else if(key.EqualsCaseInsensitive("printFreq") )	   
	    printFreq = value.ConvertToInt();
	  else if (key.EqualsCaseInsensitive("asdsfIgnoreFreq"))
	    asdsfIgnoreFreq = value.ConvertToDouble(); 
	  else if (key.EqualsCaseInsensitive("asdsfConvergence"))
	    asdsfConvergence = value.ConvertToDouble();
	  else if (key.EqualsCaseInsensitive("heatFactor"))
	    heatFactor = value.ConvertToDouble();
	  else if(key.EqualsCaseInsensitive("swapInterval"))
	    swapInterval = value.ConvertToInt();
	  else if(key.EqualsCaseInsensitive("tuneHeat"))
	    tuneHeat = convertToBool(value);
	  else if(key.EqualsCaseInsensitive("tuneFreq"))
	    tuneFreq = value.ConvertToInt();
	  else if(key.EqualsCaseInsensitive("burninGen"))
	    burninGen = value.ConvertToInt();
	  else if(key.EqualsCaseInsensitive("burninProportion"))
	    burninProportion = value.ConvertToDouble();
	  else 	      
	    cerr << "WARNING: ignoring unknown value >"  << key << "< and >" << value <<  "<" << endl; 
	}
    }
}
