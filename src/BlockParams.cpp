#include "BlockParams.hpp"
#include "axml.h"


void BlockParams::initialize(const TreeAln &traln)
{
  NCL_BLOCKTYPE_ATTR_NAME = "PARAMS"; 

  hasDna = false; 
  hasAA = false; 

  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition =  traln.getPartition(i);
      hasDna |= ( partition->dataType == DNA_DATA); 
      hasAA |= ( partition->dataType == AA_DATA ); 
      
      assert(partition->dataType == DNA_DATA || partition->dataType == AA_DATA) ; 
    }
} 

void BlockParams::parseScheme(NxsToken& token, category_t cat, nat &idCtr)
{
  vector<bool> partAppeared(false); 

  token.GetNextToken();
  auto str = token.GetToken(false); 
  assert(str.compare("=") == 0); 
  token.GetNextToken();
  str = token.GetToken(false); 
  assert(str.compare("(") == 0); 
  token.GetNextToken();
  str = token.GetToken(false); 

  bool newLink = true; 
  RandomVariable *curVar = nullptr; 
  while(str.compare(")") != 0 )
    {
      int part = str.ConvertToInt() ; 
      if(newLink)
	{
	  RandomVariable r(cat, idCtr) ; 
	  curVar = &r; 
	  newLink = false; 
	}

      curVar->addPartition(part); 
      partAppeared[part] = true; 
      
      token.GetNextToken();
      str = token.GetToken(false);
      
      if(str.compare("/") == 0)
	{
	  parameters.push_back(*curVar); 
	  newLink = true; 
	  token.GetNextToken();
	  str = token.GetToken(false); 
	}
      else if(str.compare(",") == 0)
	{
	  token.GetNextToken();
	  str = token.GetToken(false); 
	}
      else 
	{
	  assert(str.compare(")") == 0); 
	  parameters.push_back(*curVar); 
	}
    }  

  for(int i = 0; i < numPart; ++i)
    {
      if(not partAppeared[i])
	{
	  RandomVariable r(cat, idCtr); 
	  parameters.push_back(r); 
	}
    }
}


void BlockParams::Read(NxsToken &token)
{
  DemandEndSemicolon(token, "PARAMS");
  nat idCtr = 0; 

  while(true)
    {
      token.GetNextToken();
      NxsBlock::NxsCommandResult res = HandleBasicBlockCommands(token); 

      if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
	return;
      if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND))
	{
	  NxsString str = token.GetToken(false); 

	  if(str.EqualsCaseInsensitive("stateFreq"))
	    {
	      parseScheme(token,  FREQUENCIES, idCtr); 
	    }
	  else if(str.EqualsCaseInsensitive("rateHet"))
	    {
	      parseScheme(token, RATE_HETEROGENEITY, idCtr); 
	    }
	  else if(str.EqualsCaseInsensitive("revMat"))
	    {
	      parseScheme(token, SUBSTITUTION_RATES, idCtr);
	    }
	  else if(str.EqualsCaseInsensitive("branchLength"))
	    {
	      assert(NOT_IMPLEMENTED);
	    }
	  else 
	    {
	      cerr << "parsing error at " << str  << endl; 
	      exit(0);
	    }
	}
    }
}
