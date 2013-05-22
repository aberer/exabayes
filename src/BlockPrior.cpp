#include "BlockPrior.hpp"

shared_ptr<AbstractPrior> BlockPrior::parsePrior(NxsToken &token)  
{
  auto value = token.GetToken(false); 
  // cout << value << endl; 
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

  map<string,category_t> categoryNameMap = 
    {
      {"topoPr", TOPOLOGY}, 
      {"brlenpr", BRANCH_LENGTHS}, 
      {"stateFreqPr" , FREQUENCIES}, 
      {"revMatPr", SUBSTITUTION_RATES}, 
      {"shapePr", RATE_HETEROGENEITY}
    } ; 

  while(true)
    {
      token.GetNextToken();
      NxsBlock::NxsCommandResult res = HandleBasicBlockCommands(token); 

      if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
	return;
      if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND))
	{	  
	  token.GetNextToken();	  

	  auto str = token.GetToken (false); 
	  category_t cat = categoryNameMap[str]; 


	  int priorPartition = -1;  
	  if(str.compare("{") == 0)
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
	      assert(generalPriors[cat] == nullptr); 
	      generalPriors[cat] = prior; 
	    }	    
	  else 
	    {
	      map<nat, shared_ptr<AbstractPrior>>& priorsForPartition = specificPriors[cat]; 
	      assert(priorsForPartition[priorPartition] == nullptr); 
	      priorsForPartition[priorPartition] = prior; 
	    }
	}
    }  
} 
