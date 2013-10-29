#include "BlockPrior.hpp"

#include <sstream>
#include <limits>

#include "priors/DiscreteModelPrior.hpp"
#include "priors/UniformPrior.hpp"
#include "priors/ExponentialPrior.hpp"
#include "priors/DirichletPrior.hpp"
#include "priors/FixedPrior.hpp"

#include "ParallelSetup.hpp"


static void expectString( std::string expectation, NxsToken& token)
{
  bool okay = token.GetToken().EqualsCaseInsensitive(expectation.c_str()); 
  if(not okay)
    {
      std::cerr << "error while parsing the config file: expected " << expectation << " but got " << token.GetToken() << std::endl; 
      exit(-1); 
    }
}



static std::vector<double> parseValues(NxsToken &token)
{
  auto result = std::vector<double>{}; 
  
  // assumption: we have already seen the ')'

  while(token.GetToken().compare(")") != 0)
    {
      auto &&iss = std::istringstream{token.GetToken()}; 
      auto value = double{0.};
      iss >> value; 
      result.push_back(value);
      token.GetNextToken();
      if(token.GetToken().compare(",") == 0)
	token.GetNextToken();
    }
  
  return result; 
}



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
  else if(value.EqualsCaseInsensitive("disc"))
    {
      // token.GetNextToken(); 
      expectString("(", token);
      
      auto modelsProbs = std::unordered_map<ProtModel, double>{};

      auto remainder = std::numeric_limits<double>::infinity();

      while(token.GetToken().compare(")") != 0)
	{
	  token.GetNextToken();

	  auto foundRemainder = token.GetToken().EqualsCaseInsensitive("remainder"); 
	  if(foundRemainder)	// keep track of the weight 
	    {
	      if(remainder != std::numeric_limits<double>::infinity())
		{
		  std::cerr << "Encountered 'remainder' twice while defining aaPr" << std::endl; 
		  exit(-1);
		}

	      token.GetNextToken();
	      expectString("=", token); 
	      token.GetNextToken();

	      auto &&iss = std::istringstream{token.GetToken()};
	      iss >> remainder ; 
	      token.GetNextToken();
	    } 
	  else 			// simply parse that model 
	    {
	      auto modelRes = ProtModelFun::getModelFromStringIfPossible(token.GetToken()); 
	      if(not std::get<0>(modelRes) )
		{
		  std::cerr << "Error: expected " << token.GetToken() << "to be a valid protein model name" << std::endl; 
		  exit(-1); 
		}
	      auto model = std::get<1>(modelRes);

	      token.GetNextToken();
	      expectString("=", token); 

	      token.GetNextToken(); 
	      auto &&iss = std::istringstream{token.GetToken()}; 
	      auto weight = double{0.}; 
	      iss >> weight; 

	      if(modelsProbs.find(model) != modelsProbs.end())
		{
		  std::cerr << "Error: model " <<  model << "occurred more than once in your specification of a discrete amino acid model prior."  << std::endl; 
		  exit(-1); 
		}

	      modelsProbs[model] = weight; 

	      token.GetNextToken();
	    }
	}

      if(remainder != std::numeric_limits<double>::infinity())
	{
	  for(auto model : ProtModelFun::getAllModels())
	    {
	      if( modelsProbs.find(model) == modelsProbs.end()  )
		modelsProbs[model] = remainder; 
	    }
	}

      return shared_ptr<AbstractPrior>(new DiscreteModelPrior(modelsProbs));
    }
  else if(value.EqualsCaseInsensitive("dirichlet"))
    {      
      auto alphas = parseValues(token); 
      return shared_ptr<AbstractPrior>(new DirichletPrior(alphas)); 
    }
  else if(value.EqualsCaseInsensitive("fixed"))
    { 
      token.GetNextToken();
      if(token.GetToken().EqualsCaseInsensitive("empirical"))
	{
	  cerr << "not implemented yet " << endl; // TODO 
	  ParallelSetup::genericExit(-1); 
	  return nullptr; 
	} 
      else
	{
	  auto res = ProtModelFun::getModelFromStringIfPossible(token.GetToken()); 
	  auto foundProt = std::get<0>(res); 

	  if( foundProt )
	    {
	      auto model = std::get<1>(res);

	      token.GetNextToken();
	      assert(token.GetToken().compare(")" ) == 0 ); 

	      auto result = std::shared_ptr<AbstractPrior>(new DiscreteModelPrior( { {model, 1.}  } ));
	      return result;
	    }
	  else 
	    {
	      // auto fixedValues = std::vector<double>{}; 
	      auto fixedValues = parseValues(token);

	      std::cout << "found fixed values " << fixedValues << std::endl; 

	      return shared_ptr<AbstractPrior>(new FixedPrior(fixedValues));
	    }
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
      ParallelSetup::genericExit(-1); 
    }

  return nullptr; 
}


void BlockPrior::Read(NxsToken &token) 
{
  DemandEndSemicolon(token, "PRIOR");

  while(true)
    {
      token.GetNextToken();
      auto  res = HandleBasicBlockCommands(token); 

      if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
	return;
      if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND))
	{
	  auto str = token.GetToken(false).ToUpper(); 
	  auto cat = CategoryFuns::getCategoryByPriorName(str); 
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
	      std::unordered_map<nat,std::shared_ptr<AbstractPrior> > &priorsForPartition = specificPriors[int(cat)]; // BAD
	      assert(priorsForPartition[priorPartition] == nullptr); 
	      priorsForPartition[priorPartition] = prior; 
	    }
	}
    }  
} 
