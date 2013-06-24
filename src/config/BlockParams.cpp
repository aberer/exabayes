#include "BlockParams.hpp"
#include "axml.h"
#include <algorithm>

void BlockParams::parseScheme(NxsToken& token, Category cat, nat &idCtr)
{
  vector<bool> partAppeared(traln->getNumberOfPartitions(), false); 

  token.GetNextToken();
  auto str = token.GetToken(false); 
  assert(str.compare("=") == 0); 
  token.GetNextToken();
  str = token.GetToken(false); 
  assert(str.compare("(") == 0); 
  token.GetNextToken();
  str = token.GetToken(false); 

  bool newLink = true; 
  RandomVariablePtr curVar; 
  while(str.compare(")") != 0 )
    {
      int part = str.ConvertToInt() ; 
      if(newLink)
	{
	  curVar = RandomVariablePtr(new RandomVariable(cat,idCtr)); 
	  newLink = false; 
	}

      curVar->addPartition(part); 
      partAppeared[part] = true; 
      
      token.GetNextToken();
      str = token.GetToken(false);
      
      if(str.compare("/") == 0)
	{
	  parameters.push_back(curVar); 
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
	  parameters.push_back(curVar); 
	}
    }  
  
  // maybe this should not be done here /= 
  // too clumsy -- must be tested 
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      if(not partAppeared[i] && 
	 (   ( cat == Category::AA_MODEL && traln->getPartition(i)->dataType == AA_DATA ) 
	     || (cat == Category::SUBSTITUTION_RATES && traln->getPartition(i)->dataType == DNA_DATA) 
	     || ( cat != Category::SUBSTITUTION_RATES && cat != Category::AA_MODEL) ))
	{
	  
	  RandomVariablePtr r(new RandomVariable(cat, idCtr)); 
	  r->addPartition(i);
	  parameters.push_back(r); 
	}
    }
}


void BlockParams::Read(NxsToken &token)
{
  DemandEndSemicolon(token, "PARAMS");
  nat idCtr = 0; 

  std::set<Category> catsFound; 

  while(true)
    {
      token.GetNextToken();
      NxsBlock::NxsCommandResult res = HandleBasicBlockCommands(token); 

      if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
	return;

      if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND))
	{	  
	  NxsString str = token.GetToken(false); 

	  Category cat = CategoryFuns::getCategoryFromLinkLabel(str); 	  
	  parseScheme(token, cat, idCtr); 

	  
	  if(catsFound.find(cat) != catsFound.end())
	    {
	      cerr << "parsing error: found a linking scheme for category  " <<  CategoryFuns::getLongName(cat) << " twice. Aborting." ; 
	      exit(1); 
	    }

	  if(cat == Category::RATE_HETEROGENEITY)
	    {
	      string copy =  str; 
	      std::transform(copy.begin(), copy.end(), copy.begin(), ::tolower); 
	    }

	  if(cat == Category::AA_MODEL
	     || cat == Category::BRANCH_LENGTHS
	     || cat == Category::TOPOLOGY)
	    {
	      cerr <<  "not implemented"; 
	      assert(0); 
	    }	    
	}
    }
}
