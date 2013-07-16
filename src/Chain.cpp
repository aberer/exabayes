#include <sstream>
#include <unordered_map>

#include "Chain.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "ProposalFunctions.hpp"
#include "LikelihoodEvaluator.hpp" 


Chain:: Chain(randKey_t seed, std::shared_ptr<TreeAln> _traln, const std::vector<std::unique_ptr<AbstractProposal> > &_proposals, std::shared_ptr<LikelihoodEvaluator> eval) 
  : traln(_traln)
  , deltaT(0)
  , runid(0)
  , tuneFrequency(100)
  , hastings(0)
  , currentGeneration(0)
  , couplingId(0)
  , chainRand(seed)
  , relWeightSum(0)
  , bestState(std::numeric_limits<double>::lowest())
  , evaluator(eval)
{
  for(auto &p : _proposals)
    {
      std::unique_ptr<AbstractProposal> copy(p->clone()); 
      assert(copy->getPrimVar().size() != 0); 
      proposals.push_back(std::move(copy)); 
    }

  eval->evaluateFullNoBackup(*traln);  ; 
  
  const std::vector<AbstractParameter*> vars = extractVariables(); 
  prior.initialize(*traln, vars);
  // saving the tree state 
  suspend(); 
}


Chain::Chain( Chain&& rhs)   
  : traln(rhs.traln)
  , deltaT(rhs.deltaT)
  , runid(rhs.runid)
  , tuneFrequency(rhs.tuneFrequency)
  , currentGeneration(rhs.currentGeneration)
  , couplingId(rhs.couplingId) 
  , chainRand(rhs.chainRand) 
  , bestState(rhs.bestState)
  , evaluator(rhs.evaluator)
{
  for(auto &p : rhs.proposals )
    proposals.emplace_back(std::move(p->clone())); 
  prior.initialize(*traln, extractVariables()); 
  
  suspend();
  
  // cout << *traln<< endl; 
}


Chain& Chain::operator=(Chain rhs)
{
  std::swap(*this, rhs); 
  return *this; 
}


std::ostream& Chain::addChainInfo(std::ostream &out) const 
{
  return out << "[run=" << runid << ",heat=" << couplingId << ",gen=" << getGeneration() << "]" ; 
}


/**
   @brief returns the inverse temperature for this chain
 */ 
double Chain::getChainHeat()
{
  double tmp = 1. + ( deltaT * couplingId ) ; 
  double inverseHeat = 1.f / tmp; 
  assert(couplingId == 0 || inverseHeat < 1.); 
  return inverseHeat; 
}




void Chain::resume()  
{
  // std::cout << "trying to resume" << std::endl; 
  
  auto vs = extractVariables(); 

  // topology must be dealt with first   
  sort(vs.begin(), vs.end(), [] (const AbstractParameter* a,const AbstractParameter* b ){ return a->getCategory() == Category::TOPOLOGY ;   } ); 
  for(auto &v : vs)
    v->applyParameter(*traln, savedContent[v->getId()]); 

  evaluator->evaluateFullNoBackup(*traln); 

  if(fabs(likelihood - traln->getTr()->likelihood) > ACCEPTED_LIKELIHOOD_EPS)
    {
      std::cerr << "While trying to resume chain: previous chain liklihood larger than " <<
	"evaluated likelihood. This is a programming error." << std::endl; 
      std::cerr << "prev=" << likelihood << "\tnow=" << traln->getTr()->likelihood << std::endl; 

      assert(0);       
    }  

  prior.initialize(*traln, vs);
  double prNow = prior.getLnPrior(); 
  
  if(fabs(lnPr - prNow) > ACCEPTED_LNPR_EPS)
    {
      std::cerr << "While trying to resume chain: previous log prior for chain larger than" <<
	"re-evaluated prior. This is a programming error." << std::endl; 
      assert(0);       
    }

  // prepare proposals 
  relWeightSum = 0; 
  for(auto& elem : proposals)
    relWeightSum +=  elem->getRelativeWeight();
}


void Chain::printProposalState(std::ostream& out ) const 
{
  std::map<Category, std::vector<AbstractProposal*> > sortedProposals; 
  for(auto& p : proposals)
    sortedProposals[p->getCategory()].push_back(p.get()) ; 
  
  for(auto &n : CategoryFuns::getAllCategories())
    {       
      Category cat = n; 

      bool isThere = false; 
      for(auto &p : proposals)
	isThere |= p->getCategory() == cat ; 
      
      if(isThere)
	{
	  tout << CategoryFuns::getLongName(n) << ":\t";
	  for(auto &p : proposals)
	    {
	      if(p->getCategory() == cat)
		{
		  p->printNamePartitions(tout); 
		  tout  << ":"  << p->getSCtr() << "\t" ; 	      
		}
	    }
	  tout << std::endl; 
	}
    }
}


AbstractProposal* Chain::drawProposalFunction()
{ 
  double r = relWeightSum * chainRand.drawRandDouble01();   

  for(auto& c : proposals)
    {
      double w = c->getRelativeWeight(); 
      if(r < w )
	return c.get();
      else 
	r -= w; 
    }

  assert(0); 
  return proposals[0].get(); 
}


void Chain::switchState(Chain &rhs)
{
  std::swap(couplingId, rhs.couplingId); 
  std::swap(bestState, rhs.bestState); 
  // swap(chainRand, rhs.chainRand); 
  std::swap(proposals, rhs.proposals); 
}


void Chain::step()
{
#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*traln, extractVariables());
#endif

  currentGeneration++; 
  tree *tr = traln->getTr();   
  evaluator->imprint(*traln);

  // inform the rng that we produce random numbers for generation x  
  chainRand.rebase(currentGeneration);

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = getChainHeat();

  auto pfun = drawProposalFunction();
 
  /* reset proposal ratio  */
  hastings = 0; 

#ifdef DEBUG_ACCEPTANCE
  double oldPrior = prior.getLnPrior();
#endif

  // tout << pfun->getName() << std::endl; 

  pfun->applyToState(*traln, prior, hastings, chainRand);
  pfun->evaluateProposal(*evaluator, *traln, prior);

  // tout << hastings << std::endl; 
  
  double priorRatio = prior.getLnPriorRatio();
  double lnlRatio = tr->likelihood - prevLnl; 

  double testr = chainRand.drawRandDouble01();
  double acceptance = exp(( priorRatio + lnlRatio) * myHeat + hastings) ; 

#ifdef DEBUG_ACCEPTANCE
  tout  << endl << "(" << oldPrior<<  " - " << prior.getLnPriorRatio()  << ") + ( " << tr->likelihood << " - "<< prevLnl   << ") * " << myHeat << " + " << hastings << endl; 
#endif

  bool wasAccepted  = testr < acceptance; 
  
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  {
  double lnl = tr->likelihood; 
  addChainInfo(tout); 
  tout << "\t" << (wasAccepted ? "ACC" : "rej" )  << "\t"<< pfun->getName() << "\t" << std::setprecision(2) << std::fixed << lnl << "\tdelta=" << lnlRatio << "\t" << hastings << std::endl; 
  }
#endif
  // debug_printAccRejc( pfun, wasAccepted, tr->likelihood, prior.getLnPrior(), hastings); 

  if(wasAccepted)
    {
      pfun->accept();      
      prior.accept();
      if(bestState < traln->getTr()->likelihood  )
	bestState = traln->getTr()->likelihood; 
      likelihood = traln->getTr()->likelihood; 
    }
  else
    {
      pfun->resetState(*traln, prior);
      pfun->reject();
      prior.reject();

      evaluator->resetToImprinted(*traln);
    }

#ifdef DEBUG_LNL_VERIFY
  evaluator->expensiveVerify(*traln); 
#endif
  
#ifdef DEBUG_TREE_LENGTH  
  assert( fabs (traln->getTreeLengthExpensive() - traln->getTreeLength())  < 1e-6); 
#endif

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*traln, extractVariables());
#endif

  if(this->tuneFrequency <  pfun->getNumCallSinceTuning() ) 
    pfun->autotune();
}


void Chain::suspend()  
{
  auto variables = extractVariables();
  nat maxV = 0; 
  for(auto &v : variables)
    if(maxV < v->getId()) 
      maxV = v->getId(); 
  // tout << "highest id: " << maxV << std::endl; 
  savedContent.resize(maxV + 1 ); 
  for(auto& v : variables)
    savedContent[v->getId()] = v->extractParameter(*traln);

  likelihood = traln->getTr()->likelihood; 
  lnPr = prior.getLnPrior();
  
  // addChainInfo(tout);
  // tout << " SUSPEND " << likelihood << std::endl; 
}


// TODO not too much thought went into this  
namespace std 
{
  template<> class hash<AbstractParameter*>
  {
  public: 
    size_t operator()(const AbstractParameter* rhs) const 
    {
       return rhs->getId(); 
    }
  } ;

  template<>
  struct equal_to<AbstractParameter*>
  {
    bool operator()(const AbstractParameter* a, const AbstractParameter* b ) const 
    {
      return a->getId() == b->getId(); 
    }
  }; 
}




const std::vector<AbstractParameter*> Chain::extractVariables() const 
{
  std::unordered_set<AbstractParameter*> result; 
  for(auto &p : proposals)
    {
      for(auto &v : p->getPrimVar())
	result.insert(v); 
      for(auto &v : p->getSecVar())
	result.insert(v); 
    }
  
  std::vector<AbstractParameter*> result2 ; 
  result2.resize(result.size()); 
  
  for(auto &v : result)      
    result2[v->getId()] = v ; 

  return result2; 
}


const std::vector<AbstractProposal*> Chain::getProposalView() const 
{
  std::vector<AbstractProposal*> result;  
  for(auto &elem: proposals)
    result.push_back(elem.get()); 
  return result; 
}


std::ostream& operator<<(std::ostream& out, const Chain &rhs)
{
  rhs.addChainInfo(out); 
  out << "\tLnL: " << rhs.getLnLikelihood() << "\tLnPr: " << rhs.getLnPrior(); 
  return out; 
}


void Chain::sample( const TopologyFile &tFile, const ParameterFile &pFile  ) const
{
  tFile.sample( *traln, getGeneration() ); 
  pFile.sample( *traln, extractVariables(), getGeneration(), prior.getLnPrior()); 
}


void Chain::readFromCheckpoint( std::ifstream &in ) 
{
  chainRand.readFromCheckpoint(in); 
  in >> couplingId; 
  readDelimiter(in);
  in >> likelihood; 
  readDelimiter(in); 
  in >> lnPr; 
  readDelimiter(in); 
  in >> currentGeneration; 
  readDelimiter(in); 
  
  std::unordered_map<std::string, AbstractProposal*> name2proposal; 
  for(auto &p :proposals)
    {
      std::stringstream ss; 
      p->printShort(ss); 
      assert(name2proposal.find(ss.str()) == name2proposal.end()); // not yet there 
      name2proposal[ss.str()] = p.get();
    }

  nat ctr = 0; 
  while(ctr < proposals.size())
    {
      std::string name; 
      std::getline(in, name, DELIM); 
      if(name2proposal.find(name) == name2proposal.end())
	{
	  std::cerr << "Could not parse the checkpoint file.  A reason for this may be that\n"
		    << "you used a different configuration or alignment file in combination\n"
		    << "with this checkpoint file. Fatality." << std::endl; 
	  exit(0); 
	}
      name2proposal[name]->readFromCheckpoint(in);
      ++ctr; 
    }


  std::unordered_map<std::string, AbstractParameter*> name2parameter; 
  for(auto &p : extractVariables())
    {
      std::stringstream ss;       
      p->printShort(ss);
      assert(name2parameter.find(ss.str()) == name2parameter.end()); // not yet there 
      name2parameter[ss.str()] = p;
    }  

  ctr = 0; 
  while(ctr < name2parameter.size())
    {
      std::string name; 
      std::getline(in, name, DELIM); 
      if(name2parameter.find(name) == name2parameter.end())
	{
	  std::cerr << "Could not parse the checkpoint file. A reason for this may be that\n"
		    << "you used a different configuration or alignment file in combination\n"
		    << "with this checkpoint file. Fatality." << std::endl; 
	  exit(0); 
	}
      auto param  = name2parameter[name]; 
      auto content = param->extractParameter(*traln); // initializes the object correctly. the object must "know" how many values are to be extracted 
      content.readFromCheckpoint(in);
      param->applyParameter(*traln, content); 
      ++ctr;
    }  
  
  // resume the chain and thereby assert that everything could be
  // restored correctly
  auto vars = extractVariables(); 
  savedContent.resize(vars.size()); 
  for(auto &v : vars)
    savedContent[v->getId()] = v->extractParameter(*traln);
  resume();
}
 
void Chain::writeToCheckpoint( std::ofstream &out) const
{
  chainRand.writeToCheckpoint(out);   
  out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 2 );
  out << couplingId << DELIM; 
  out << likelihood << DELIM; 
  out << lnPr << DELIM; 
  out << currentGeneration << DELIM; // well ... 

  for(auto &p : proposals)
    p->writeToCheckpoint(out);
  
  for(auto &var: extractVariables())
    {
      auto compo = var->extractParameter(*traln); 
      var->printShort(out); 
      out << DELIM; 
      compo.writeToCheckpoint(out);
    }  
}   
