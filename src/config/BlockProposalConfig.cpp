#include <cassert>

#include "BlockProposalConfig.hpp"
#include "comm/ParallelSetup.hpp"


extern void genericExit(int code); 


BlockProposalConfig::BlockProposalConfig()
  : etbrStopProb(0.5)
  , esprStopProp(0.5)    
  , parsimonyWarp(0.1)
  , parsSPRRadius(-1)
{
  NCL_BLOCKTYPE_ATTR_NAME = "PROPOSALS"; 
}


void BlockProposalConfig::Read(NxsToken &token)
{   
  DemandEndSemicolon(token, "PROPOSALS");
  
  auto ps = ProposalTypeFunc::getAllProposals(); 

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

	  if(ProposalTypeFunc::isValidName(key))
	    {
	      double val = value.ConvertToDouble(); 
	      auto t = ProposalTypeFunc::getTypeFromConfigString(key); 
	      if(userValue.find(t) != userValue.end())
		{
		  std::cerr << "encountered the value " << key << "twice in the config file" << std::endl; 
		  ParallelSetup::genericExit(-1); 
		}
	      else 
		userValue[t] = val; 
	    }
	  else if(key.EqualsCaseInsensitive("esprstopprob"))	    
	    esprStopProp = value.ConvertToDouble();	  
	  else if(key.EqualsCaseInsensitive("etbrstopprob"))
	    etbrStopProb = value.ConvertToDouble(); 	  
	  else if(key.EqualsCaseInsensitive("guidedsprradius"))
	    guidedRadius = value.ConvertToInt();
	  else if(key.EqualsCaseInsensitive("parsimonyWarp"))	    
	    parsimonyWarp = value.ConvertToDouble();
	  else if(key.EqualsCaseInsensitive("parssprradius"))
	    {
	      parsSPRRadius = value.ConvertToInt();
	      // tout << "\n\nfound spr radius " << parsSPRRadius << "\n\n" << std::endl ;
	    }
	  else 	      
	    {
	      cerr << "WARNING: ignoring unknown value >"  << key << "< and >" << value <<  "<" << endl; 
	      assert(0);
	    }
	}
    }
}


void BlockProposalConfig::verify()
{
  if(not (0.01 < esprStopProp && esprStopProp <= 1 ))
    {
      tout << "Error: >esprStopProp< must be in the range (0.01,1]" << std::endl; 
      ParallelSetup::genericExit(-1); 
    }

  if(not (0.01 < esprStopProp && esprStopProp <= 1 ))
    {
      tout << "Error: >etbrStopProb< must be in the range (0.01,1]" << std::endl; 
      ParallelSetup::genericExit(-1); 
    }
  
  if(not (0.001 < parsimonyWarp && parsimonyWarp < 10))
    {
      tout << "Error: >parsimonyWarp< must be in the range (0.001,10)" << std::endl; 
      ParallelSetup::genericExit(-1); 
    }
  
  if(parsSPRRadius != -1 &&  (parsSPRRadius <= 1 ))
    {
      tout << "Error: >parsSPRRadius< must be in the range (1,inf)" << std::endl; 
      ParallelSetup::genericExit(-1); 
    }
}

// NOTICE 

