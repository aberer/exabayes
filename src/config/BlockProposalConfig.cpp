#include <cassert>

#include "BlockProposalConfig.hpp"

BlockProposalConfig::BlockProposalConfig()
{
  NCL_BLOCKTYPE_ATTR_NAME = "PROPOSALS"; 
  // setupMap();
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
		  exit(0); 
		}
	      else 
		userValue[t] = val; 
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

