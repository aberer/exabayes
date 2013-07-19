#include "BlockParams.hpp"
#include "axml.h"
#include "GlobalVariables.hpp"
#include "ParallelSetup.hpp"

#include <algorithm>

extern void genericExit(int code); 


void BlockParams::partitionError(nat partition, nat totalPart) const
{
  std::cerr << "In the parameter block of the configuration file you specified partition " << partition << ". However, there are only " << totalPart << " partitions in total in your alignment." << std::endl; 
  ParallelSetup::genericExit(-1); 
}


void BlockParams::parseScheme(NxsToken& token, Category cat, nat &idCtr)
{
  nat numPart = traln->getNumberOfPartitions();
  vector<bool> partAppeared(numPart, false); 

  token.GetNextToken();
  auto str = token.GetToken(false); 
  assert(str.compare("=") == 0); 
  token.GetNextToken();
  str = token.GetToken(false); 
  assert(str.compare("(") == 0); 
  token.GetNextToken();
  str = token.GetToken(false); 

  bool newLink = true; 
  std::unique_ptr<AbstractParameter> curVar; 
  while(str.compare(")") != 0 )
    {
      nat part = str.ConvertToInt() ; 
      if(newLink)
	{
	  curVar = CategoryFuns::getParameterFromCategory(cat, idCtr);
	  idCtr++; 
	  newLink = false; 
	}

      if(not (part < numPart ))
	partitionError(part, numPart); 

      curVar->addPartition(part); 
      partAppeared[part] = true; 
      
      token.GetNextToken();
      str = token.GetToken(false);
      
      if(str.compare("/") == 0)
	{
	  parameters.push_back(std::move(curVar)); 
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
	  parameters.push_back(std::move(curVar)); 
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
	  auto r = CategoryFuns::getParameterFromCategory(cat, idCtr);
	  ++idCtr; 

	  r->addPartition(i);
	  parameters.push_back(std::move(r)); 
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
	      ParallelSetup::genericExit(-1); 
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


vector<unique_ptr<AbstractParameter> > BlockParams::getParameters() const
{
  vector<unique_ptr<AbstractParameter> > result; 
  for(auto &p : parameters	)
    result.push_back(std::unique_ptr<AbstractParameter>(p->clone() )); 
  return result; 
}
