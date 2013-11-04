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
  nat numPart = tralnPtr->getNumberOfPartitions();
  auto partAppeared = std::vector<bool>(numPart, false); 

  token.GetNextToken();
  assert(token.GetToken().EqualsCaseInsensitive("=")); 
  token.GetNextToken();
  assert(token.GetToken().EqualsCaseInsensitive("(")); 

  auto scheme = std::vector<std::vector<nat>>{}; 
  while(not token.GetToken().EqualsCaseInsensitive(")") )
    {
      auto schemePart = std::vector<nat>{}; 
      bool startingNext = true; 
      
      // parse one partition part 
      while( not token.GetToken().EqualsCaseInsensitive(")") 	     
	     && (startingNext || not token.GetToken().EqualsCaseInsensitive(",")) )
	{
	  startingNext = false; 
	  // check if the separation item is an expasion item (-)
	  auto doExpand = token.GetToken().EqualsCaseInsensitive("-"); 
	  token.GetNextToken();
	  nat part = token.GetToken().ConvertToInt();
	  nat start = doExpand ? schemePart.back() +1  : part; 
	  nat end = part + 1 ;  
	  for(nat i = start ;  i < end ; ++i )
	    {
	      if(not (i < numPart))
		partitionError(i, numPart); 
	      if(partAppeared.at(i))
		{
		  tout << "error: partition " << i << " occurring twice in the same scheme. Check your parameter-block!" << std::endl; 
		  exit(-1); 
		}
	      partAppeared.at(i) = true; 
	      schemePart.push_back(i); 
	    }
	  
	  token.GetNextToken(); 
	}

      scheme.push_back(schemePart); 
    }

  // instantiate the parameters 
  for(auto schemePart : scheme )
    {
      assert(schemePart.size()  > 0 ); 
      parameters.push_back(CategoryFuns::getParameterFromCategory(cat,idCtr,getNumSeen(cat), schemePart));
      ++idCtr; 
      tout << "initialized " << parameters.back().get() << std::endl; 
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
	  auto str = token.GetToken(false); 

	  Category cat = CategoryFuns::getCategoryFromLinkLabel(str); 	  
	  parseScheme(token, cat, idCtr); 

	  if(catsFound.find(cat) != catsFound.end())
	    {
	      cerr << "parsing error: found a linking scheme for category  " <<  CategoryFuns::getLongName(cat) << " twice. Aborting." ; 
	      ParallelSetup::genericExit(-1); 
	    }

	  if( cat == Category::TOPOLOGY)
	    {
	      cerr <<  "not implemented"; 
	      assert(0); 
	    }	    
	}
    }
}


std::vector<std::unique_ptr<AbstractParameter> > BlockParams::getParameters() const
{
  std::vector<std::unique_ptr<AbstractParameter> > result; 
  for(auto &p : parameters )
    result.push_back(std::unique_ptr<AbstractParameter>(p->clone() )); 
  return result; 
}



nat BlockParams::getNumSeen(Category cat) 
{
  nat result = 0; 
  for(auto &p : parameters)
    if(p->getCategory() == cat )
      ++result;
  return result; 
} 
