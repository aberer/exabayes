#include "BlockPrior.hpp"


shared_ptr<AbstractPrior> BlockPrior::parsePrior(NxsToken &token)  
{
  auto value = token.GetToken(false); 
  token.GetNextToken();

  assert(token.GetToken(false).compare("(") == 0); 

  if(value.EqualsCaseInsensitive("uniform"))
    {
      token.GetNextToken();

      // for non-continuous variables (e.g., topology)
      if(token.GetToken().compare(")") == 0) 
	return shared_ptr<AbstractPrior> (new UniformPrior(0,0)); 

      double n1 = atof(token.GetToken().c_str()); 
      token.GetNextToken();
      assert(token.GetToken().compare(",") == 0);
      token.GetNextToken();
      double n2 = atof(token.GetToken().c_str());
      token.GetNextToken();
      assert(token.GetToken().compare(")") == 0);
      return shared_ptr<AbstractPrior> (new UniformPrior(n1,n2));  
    }
  else if(value.EqualsCaseInsensitive("dirichlet"))
    {      
      vector<double> alphas; 
      while(token.GetToken().compare(")") != 0 )
	{
	  token.GetNextToken();
	  alphas.push_back(atof(token.GetToken().c_str())); 
	  token.GetNextToken();
	  assert(token.GetToken().compare(",") == 0 
		 || token.GetToken().compare(")") == 0); 
	} 
      return shared_ptr<AbstractPrior>(new DirichletPrior(alphas)); 
    }
  else if(value.EqualsCaseInsensitive("fixed"))
    { 
      token.GetNextToken();
      if(token.GetToken().EqualsCaseInsensitive("empirical"))
	{
	  cerr << "not implemented yet " << endl; // TODO 
	  exit(1);
	  return nullptr; 
	} 
      else 
	{
	  vector<double> fixedValues; 	  
	  while(token.GetToken().compare(")") != 0)
	    {
	      fixedValues.push_back(atof(token.GetToken().c_str()));
	      token.GetNextToken();
	      if(token.GetToken().compare(",") == 0)
		token.GetNextToken();
	    }
	  return shared_ptr<AbstractPrior>(new FixedPrior(fixedValues));
	}
    }
  else if(value.EqualsCaseInsensitive("exponential"))
    {
      token.GetNextToken();
      double n1 = atof(token.GetToken().c_str());
      token.GetNextToken();
      assert(token.GetToken().compare(")") == 0);
      return shared_ptr<AbstractPrior>(new ExponentialPrior(n1));
    }
  else 
    {
      cerr << "attempted to parse prior. Did not recognize keyword " <<  value << endl; 
      exit(1);
    }

  return nullptr; 
}


void BlockPrior::Read(NxsToken &token) 
{
  DemandEndSemicolon(token, "PRIOR");

  while(true)
    {
      token.GetNextToken();
      NxsBlock::NxsCommandResult res = HandleBasicBlockCommands(token); 

      if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
	return;
      if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND))
	{
	  auto str = token.GetToken(false).ToUpper(); 
	  Category cat = CategoryFuns::getCategoryByPriorName(str); 
	  token.GetNextToken();

	  int priorPartition = -1;  
	  if(token.GetToken().compare("{") == 0)
	    {
	      token.GetNextToken();
	      str = token.GetToken(false);
	      priorPartition = str.ConvertToInt();
	      token.GetNextToken();
	      str = token.GetToken(false); 
	      assert(str.compare("}") == 0) ;
	      token.GetNextToken();
	    }

	  assert(priorPartition < (int)numPart); 

	  auto prior = parsePrior(token);

	  if(priorPartition == -1 )	  
	    {
	      assert(generalPriors[int(cat)] == nullptr); // BAD
	      generalPriors[int(cat)] = prior;		  // BAD
	    }	    
	  else 
	    {
	      map<nat,PriorPtr> &priorsForPartition = specificPriors[int(cat)]; // BAD
	      assert(priorsForPartition[priorPartition] == nullptr); 
	      priorsForPartition[priorPartition] = prior; 
	    }
	}
    }  
} 
