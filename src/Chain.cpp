#include <sstream>

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


Chain::Chain( const Chain& rhs)   
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



std::ostream& Chain::addChainInfo(std::ostream &out) const 
{
  return out << "[run " << runid << ",heat " << couplingId << "]" ; 
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
  auto vs = extractVariables(); 

  // topology must be dealt with first   
  sort(vs.begin(), vs.end(), [] (const AbstractParameter* a,const AbstractParameter* b ){ return a->getCategory() == Category::TOPOLOGY ;   } ); 
  for(auto & v : vs)
    v->applyParameter(*traln, v->getSavedContent()); 

  evaluator->evaluateFullNoBackup(*traln); 

  if(fabs(likelihood - traln->getTr()->likelihood) > ACCEPTED_LIKELIHOOD_EPS)
    {
      std::cerr << "While trying to resume chain: previous chain liklihood larger than " <<
	"evaluated likelihood. This is a programming error." << std::endl; 
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


void Chain::debug_printAccRejc(AbstractProposal *prob, bool accepted, double lnl, double lnPr ) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  tout << "[run=" << runid << ",heat="
       << couplingId << ",gen="  << currentGeneration << "]\t" 
       << (accepted ? "ACC" : "rej" )  << "\t"<< prob->getName() << "\t"
       << std::setprecision(2) << std::fixed << lnl  << std::endl; //   
#endif
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
  debug_printAccRejc( pfun, wasAccepted, tr->likelihood, prior.getLnPrior()); 

  if(wasAccepted)
    {
      pfun->accept();      
      prior.accept();
      if(bestState < traln->getTr()->likelihood  )
	bestState = traln->getTr()->likelihood; 
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
    {
      pfun->autotune();
    }

}


void Chain::suspend()  
{
  auto variables = extractVariables();
  for(auto & v : variables)    
    v-> setSavedContent(v->extractParameter(*traln)) ; 
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
