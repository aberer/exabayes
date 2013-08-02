#include <sstream>
#include <map> 
#include <unordered_map>

#include "Chain.hpp"		
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"

#include "LikelihoodEvaluator.hpp" 
#include "ParallelSetup.hpp"
#include "Category.hpp"

// #define VERIFY_GEN 4
#define VERIFY_GEN 1000000


Chain:: Chain(randKey_t seed, std::shared_ptr<TreeAln> _traln, 
	      const std::vector<std::unique_ptr<AbstractProposal> > &_proposals, 
	      std::vector<ProposalSet> _proposalSets, 
	      std::unique_ptr<LikelihoodEvaluator> eval) 
  : tralnPtr(_traln)
  , deltaT(0)
  , runid(0)
  , tuneFrequency(100)
  , hastings(0)
  , currentGeneration(0)
  , couplingId(0)
  , proposalSets(_proposalSets)
  , chainRand(seed)
  , relWeightSumSingle(0)
  , relWeightSumSets(0)
  , bestState(std::numeric_limits<double>::lowest())
  , evaluator(std::move(eval))
{
  for(auto &p : _proposals)
    {
      std::unique_ptr<AbstractProposal> copy(p->clone()); 
      assert(copy->getPrimaryParameterView().size() != 0); 
      proposals.push_back(std::move(copy)); 
    }

  Branch root(tralnPtr->getTr()->start->number, tralnPtr->getTr()->start->back->number); 
  evaluator->evaluate(*tralnPtr, root, true); 
  
  const std::vector<AbstractParameter*> vars = extractParameters(); 
  prior.initialize(*tralnPtr, vars);

  // saving the tree state 
  suspend(); 
}


Chain::Chain( Chain&& rhs)   
  : tralnPtr(rhs.tralnPtr)
  , deltaT(rhs.deltaT)
  , runid(rhs.runid)
  , tuneFrequency(rhs.tuneFrequency)
  , currentGeneration(rhs.currentGeneration)
  , couplingId(rhs.couplingId) 
  , proposalSets(rhs.proposalSets)
  , chainRand(rhs.chainRand) 
  , bestState(rhs.bestState)
  , evaluator(std::move(rhs.evaluator))
{
  for(auto &p : rhs.proposals )
    proposals.emplace_back(std::move(p)); // ->clone
  prior.initialize(*tralnPtr, extractParameters()); 
  suspend();
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


void Chain::resume(bool evaluate, bool checkLnl) 
{    
  auto vs = extractParameters(); 

  // set the topology first 
  bool topoFound = false; 
  for(auto &v : vs)
    {
      if(v->getCategory( ) == Category::TOPOLOGY)	
	{
	  assert(not topoFound);
	  topoFound = true; 
	  // tout << "RESUME: APPLY " << v << "\t" << savedContent[v->getId()] << std::endl; 
	  v->applyParameter(*tralnPtr, savedContent[v->getId()]); 
	}
    }

  // now deal with all the other parameters 
  for(auto &v : vs)
    {
      if(v->getCategory() != Category::TOPOLOGY)
	{
	  // tout << "RESUME: APPLY " << v << "\t" << savedContent[v->getId()] << std::endl; 
	  v->applyParameter(*tralnPtr, savedContent[v->getId()]); 
	}
    }

  if(evaluate)
    {
      Branch root(tralnPtr->getTr()->start->number, tralnPtr->getTr()->start->back->number); 
      evaluator->evaluate(*tralnPtr, root, true);

      if(fabs(likelihood - tralnPtr->getTr()->likelihood) >  ACCEPTED_LIKELIHOOD_EPS )
	{
	  addChainInfo(std::cerr); 
	  std::cerr << "While trying to resume chain: previous chain liklihood"
		    <<  " could not be exactly reproduced. Please report this issue." << std::endl; 
	  std::cerr << MAX_SCI_PRECISION << 
	    "prev=" << likelihood << "\tnow=" << tralnPtr->getTr()->likelihood << std::endl; 	  
	  assert(0);       
	}  
      
      prior.initialize(*tralnPtr, vs);
      double prNow = prior.getLnPrior(); 
  
      if( fabs(lnPr - prNow) > ACCEPTED_LNPR_EPS)
	{
	  std::cerr << MAX_SCI_PRECISION ; 
	  std::cerr << "While trying to resume chain: previous log prior could not be\n"
		    << "reproduced. This is a programming error." << std::endl; 
	  std::cerr << "Prior was " << lnPr << "\t now we have " << prNow << std::endl; 
	  assert(0);       
	}
    }

  // prepare proposals 
  relWeightSumSingle = 0; 
  for(auto& elem : proposals)
    relWeightSumSingle +=  elem->getRelativeWeight();
  
  relWeightSumSets = 0; 
  for(auto &set : proposalSets)
    relWeightSumSets += set.getRelativeWeight();

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


ProposalSet& Chain::drawProposalSet()
{
  double r = relWeightSumSets * chainRand.drawRandDouble01(); 

  for(auto &c : proposalSets )
    {
      double w = c.getRelativeWeight(); 
      if(r < w )
	return c; 
      else 
	r -= w ; 
    }

  assert(0); 
  return proposalSets[0]; 
}


AbstractProposal* Chain::drawProposalFunction()
{ 
  double r = relWeightSumSingle * chainRand.drawRandDouble01();   

  for(auto& c : proposals)
    {
      double w = c->getRelativeWeight(); 
      if(r < w )
	{
	  // std::cout << "drawn " << c.get() << std::endl; 
	  return c.get();
	}
      else 
	r -= w; 
    }

  assert(0); 
  return proposals[0].get(); 
}


std::string
Chain::serializeConditionally( CommFlag commFlags)  const 
{
  std::stringstream ss; 
  
  if(commFlags & CommFlag::PrintStat)
    {
      ss.write(reinterpret_cast<const char*>(&couplingId), sizeof(couplingId)); 
      ss.write(reinterpret_cast<const char*>(&bestState), sizeof(bestState )); 
      ss.write(reinterpret_cast<const char*>(&likelihood), sizeof(likelihood)); 
      ss.write(reinterpret_cast<const char*>(&lnPr), sizeof(lnPr)); 
      ss.write(reinterpret_cast<const char*>(&currentGeneration), sizeof(currentGeneration)); 

      // std::cout << "WROTE "<< couplingId << std::endl; 
      // std::cout << "WROTE " << bestState << std::endl; 
      // std::cout << "WROTE " << likelihood << std::endl; 
      // std::cout << "WROTE " << lnPr << std::endl; 

    }

  if(commFlags & CommFlag::Proposals)
    {
      assert(proposals.size() > 0);  
      for(auto& p : proposals)
	p->writeToCheckpoint(ss); 
    }

  if(commFlags & CommFlag::Tree)
    {
      for(auto &var: extractParameters())
	{
	  const ParameterContent& compo = savedContent.at(var->getId()); 
	  
	  std::stringstream tmp; 
	  var->printShort(tmp); 	  	  
	  this->writeString(ss,tmp.str()); 
	  compo.writeToCheckpoint(ss);
	}
    }

  return ss.str(); 
}


void
Chain::deserializeConditionally(std::string str, CommFlag commFlags)
{    
  std::stringstream ss; 
  ss.str(str); 

  if(commFlags & CommFlag::PrintStat)
    {
      ss.read(reinterpret_cast<char*>(&couplingId), sizeof(couplingId) ); 
      ss.read(reinterpret_cast<char*>(&bestState), sizeof(bestState)); 
      ss.read(reinterpret_cast<char*>(&likelihood), sizeof(likelihood)); 
      ss.read(reinterpret_cast<char*>(&lnPr), sizeof(lnPr)); 
      ss.read(reinterpret_cast<char*>(&currentGeneration), sizeof(currentGeneration)); 
      std::cout << MAX_SCI_PRECISION; 

      // std::cout << "READ " << couplingId << std::endl; 
      // std::cout << "READ " << bestState << std::endl; 
      // std::cout << "READ " << likelihood << std::endl; 
      // std::cout << "READ " << lnPr << std::endl; 
    }

  if(commFlags & CommFlag::Proposals)
    initProposalsFromStream(ss);
  
  if(commFlags & CommFlag::Tree)
    {
      std::unordered_map<std::string, AbstractParameter*> name2parameter; 
      for(auto &p : extractParameters())
	{
	  std::stringstream ss;       
	  p->printShort(ss);
	  assert(name2parameter.find(ss.str()) == name2parameter.end()); // not yet there 
	  name2parameter[ss.str()] = p;
	}  

      nat ctr = 0; 
      while(ctr < name2parameter.size())
	{
	  // std::string name = cRead<std::string>(in) ; 
	  std::string name = readString(ss); 
	  if(name2parameter.find(name) == name2parameter.end())
	    {
	      std::cerr << "Could not parse the checkpoint file. A reason for this may be that\n"
			<< "you used a different configuration or alignment file in combination\n"
			<< "with this checkpoint file. Fatality." << std::endl; 
	      ParallelSetup::genericExit(-1); 
	    }
	  auto param  = name2parameter[name]; 

	  ParameterContent content = param->extractParameter(*tralnPtr); // initializes the object correctly. the object must "know" how many values are to be extracted 
	  content.readFromCheckpoint(ss);
	  savedContent[param->getId()]  = content; 

	  ++ctr;
	}
    }
}




void Chain::stepSingleProposal()
{
  // DEBUG 
  // tout << "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << currentGeneration << std::endl; 
  // {
  //   assert(tralnPtr->getNumberOfPartitions() == 1 ); 
  //   auto partition = tralnPtr->getPartition(0); 

  //   for(nat i = 0; i < tralnPtr->getNumberOfInnerNodes(); ++i)
  //     {
  // 	nat index = i + tralnPtr->getNumberOfTaxa( )+ 1 ; 
  // 	tout << "[ array " << index <<  "] " <<  partition->xVector[i] << std::endl; 
  //     }
  // }


  // evaluator->expensiveVerify(*tralnPtr);   
  
  // DEBUG 
  // tout << *tralnPtr << std::endl; 
  // for(int i = 0; i < tralnPtr->getNumberOfPartitions() ; ++i) 
  //   {
  //     auto p = tralnPtr->getPartition(i); 
  //     tout << *p << std::endl; 
  //   }
  // tout << "fracchange=" << tralnPtr->getTr()->fracchange << std::endl; 
  
  // if(currentGeneration > 138)
  //   exit(0); 

  // debugPrint = (  currentGeneration >= VERIFY_GEN   ) ; 


  double prevLnl = tralnPtr->getTr()->likelihood; 
  double myHeat = getChainHeat();
  
  auto pfun = drawProposalFunction();
 
  /* reset proposal ratio  */
  hastings = 0; 

  pfun->applyToState(*tralnPtr, prior, hastings, chainRand);

  if(debugPrint)
    tout << TreePrinter(false, true, false).printTree(*tralnPtr) << std::endl; 

  pfun->evaluateProposal(evaluator.get(), *tralnPtr, prior);
  
  double priorRatio = prior.getLnPriorRatio();
  double lnlRatio = tralnPtr->getTr()->likelihood - prevLnl; 

  double testr = chainRand.drawRandDouble01();
  double acceptance = exp(( priorRatio + lnlRatio) * myHeat + hastings) ; 

  bool wasAccepted  = testr < acceptance; 
  
#ifdef DEBUG_SHOW_EACH_PROPOSAL 
  addChainInfo(tout); 
  tout << "\t" << (wasAccepted ? "ACC" : "rej" )  << "\t"<< pfun->getName() << "\t" 
       << SOME_FIXED_PRECISION << prevLnl << "\tdelta(lnl)=" << lnlRatio << "\tdelta(lnPr)" << priorRatio << "\t" << hastings << std::endl; 
#endif

  if(wasAccepted)
    {
      pfun->accept();      
      prior.accept();
      if(bestState < tralnPtr->getTr()->likelihood  )
	bestState = tralnPtr->getTr()->likelihood; 
      likelihood = tralnPtr->getTr()->likelihood; 
      lnPr = prior.getLnPrior();
    }
  else
    {
      pfun->resetState(*tralnPtr, prior);
      pfun->reject();
      prior.reject();

      evaluator->resetToImprinted(*tralnPtr);
    }

  if(this->tuneFrequency <  pfun->getNumCallSinceTuning() ) 
    pfun->autotune();
}


void Chain::stepSetProposal()
{
  double myHeat = getChainHeat(); 
  auto pSet = drawProposalSet(); 

  std::vector<double> oldPartitionLnls = tralnPtr->getPartitionLnls(); 

  std::unordered_map<AbstractProposal*, double> p2Hastings; 
  std::unordered_map<AbstractProposal*, double> p2LnPriorRatio; 
  std::unordered_map<AbstractProposal*, double> p2OldLnl; 

  std::vector<nat> affectedPartitions; 

  for(auto &proposal : pSet.getProposalView())
    {
      double lHast = 0;       
      // double oldPr = prior.getLnPriorRatio();
      prior.reject();
      proposal->applyToState(*tralnPtr, prior, lHast, chainRand); 
      p2LnPriorRatio[proposal] = prior.getLnPriorRatio(); 

      p2Hastings[proposal] = lHast; 

      for(auto &p : proposal->getPrimaryParameterView())
	{
	  auto partitions = p->getPartitions(); 
	  affectedPartitions.insert(affectedPartitions.end(), partitions.begin(), partitions.end()); 
	  
	  double lnl = 0; 
	  for(auto &partition : partitions)
	    lnl += oldPartitionLnls[partition]; 
	  p2OldLnl[proposal] = lnl; 
	}
    }

  // no partition must occur twice!
  std::unordered_set<nat> uniquePartitions; 
  for(auto &p : affectedPartitions)
    {
      assert(uniquePartitions.find(p) == uniquePartitions.end());
      uniquePartitions.insert(p); 
    }

  evaluator->evaluatePartitions(*tralnPtr, affectedPartitions, true );
  std::vector<double> newPLnls = tralnPtr->getPartitionLnls();

  prior.reject();		// slight abuse 
  std::vector<nat> partitionsToReset; 
  nat accCtr = 0; 
  nat total = 0; 
  for(auto &proposal : pSet.getProposalView())
    {
      ++total; 
      double newLnl =0 ;  
      for(auto var : proposal->getPrimaryParameterView())
	for(auto p : var->getPartitions()) 
	  newLnl += newPLnls[p]; 
 
      double accRatio = exp(( p2LnPriorRatio[proposal] + ( newLnl - p2OldLnl[proposal]) ) * myHeat + p2Hastings[proposal]);

      if(chainRand.drawRandDouble01() < accRatio)
	{
	  ++accCtr;
	  proposal->accept();
	}
      else 
	{
	  // TODO prior more efficient 

	  proposal->resetState(*tralnPtr, prior);
	  proposal->reject();	  
	  for(auto& var: proposal->getPrimaryParameterView())
	    {
	      auto partitions = var->getPartitions();
	      partitionsToReset.insert(partitionsToReset.end(), partitions.begin(), partitions.end());
	    }
	}      
    }
  evaluator->resetSomePartitionsToImprinted(*tralnPtr, partitionsToReset); 

  newPLnls = tralnPtr->getPartitionLnls();
  for(auto &partition : partitionsToReset)
    newPLnls[partition] = oldPartitionLnls[partition]; 
  tralnPtr->setPartitionLnls(newPLnls); 
  double lnl = std::accumulate(newPLnls.begin(), newPLnls.end(), 0.);

  // meh 
  tralnPtr->getTr()->likelihood = lnl;
  likelihood = lnl; 
  lnPr = prior.getLnPrior();
  
  if(bestState < lnl)
    bestState = lnl; 

#ifdef DEBUG_SHOW_EACH_PROPOSAL
  addChainInfo(tout);
  tout << "\t" << accCtr << "/"  << total << "\t" << pSet << "\t" << likelihood << std::endl; 
#endif

  for(auto &proposal : pSet.getProposalView())
    if(this->tuneFrequency < proposal->getNumCallSinceTuning()) // meh
      proposal->autotune(); 

  // double meh 
  prior.initialize(*tralnPtr, extractParameters());
}


void Chain::step()
{
  currentGeneration++; 

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*tralnPtr, extractParameters());
#endif

  evaluator->imprint(*tralnPtr);
  // inform the rng that we produce random numbers for generation x  
  chainRand.rebase(currentGeneration);

  double sum = relWeightSumSingle + relWeightSumSets; 
  if(chainRand.drawRandDouble01() * sum < relWeightSumSingle)
    stepSingleProposal();
  else 
    stepSetProposal();

#ifdef DEBUG_LNL_VERIFY
  evaluator->expensiveVerify(*tralnPtr); 
#endif

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*tralnPtr, extractParameters());
#endif

  if( currentGeneration == VERIFY_GEN  )
    {      
      tout << "EVAL" << std::endl; 
      evaluator->evaluate(*tralnPtr, evaluator->findVirtualRoot(*tralnPtr), true); 
    }
}


 
void Chain::suspend()  
{
  auto variables = extractParameters();
  savedContent.clear(); 

  for(auto& v : variables)
    {
      assert(savedContent.find(v->getId()) == savedContent.end()); 
      savedContent[v->getId()] = v->extractParameter(*tralnPtr); 
    }

  // if(not paramsOnly)
  //   {
#ifdef EFFICIENT
  // too expensive 
  assert(0); 
#endif
  // resume(false, true ); 
  // Branch rootBranch(tralnPtr->getTr()->start->number, tralnPtr->getTr()->start->back->number);
      // evaluator->evaluate(*tralnPtr, rootBranch, true);
  likelihood = tralnPtr->getTr()->likelihood; 
  lnPr = prior.getLnPrior();
  // }
  
  // addChainInfo(tout); 
  // tout << " SUSPEND " << likelihood << "\t" << lnPr << std::endl; 
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




const std::vector<AbstractParameter*> Chain::extractParameters() const 
{
  std::unordered_set<AbstractParameter*> result; 

  // add parameters in the default proposals 
  for(auto &p : proposals)
    {
      for(auto &v : p->getPrimaryParameterView())
	result.insert(v); 
      for(auto &v : p->getSecondaryParameterView())
	result.insert(v); 
    }
  
  // add the proposals in the proposal sets 
  for(auto &p : proposalSets)
    {
      for(auto &p2 : p.getProposalView() )
	{
	  auto tmp =  p2->getPrimaryParameterView();
	  result.insert( tmp.begin(), tmp.end()); 
	  tmp = p2->getSecondaryParameterView();
	  result.insert(tmp.begin(), tmp.end()); 
	}
    }

  std::vector<AbstractParameter*> result2 ;   
  for(auto &v : result) 
    result2.push_back(v); 

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
  // out  << "\tLnl: " << rhs.tralnPtr->getTr()->likelihood << "\tLnPr: " << rhs.prior.getLnPrior();
  out  << "\tLnl: " << rhs.likelihood << "\tLnPr: " << rhs.lnPr; 
  return out;  
}


void Chain::sample(  TopologyFile &tFile,  ParameterFile &pFile  )  const
{
  tFile.sample( *tralnPtr, getGeneration() ); 
  pFile.sample( *tralnPtr, extractParameters(), getGeneration(), prior.getLnPrior()); 
}



void Chain::initProposalsFromStream(std::istream& in)
{
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
      std::string name = readString(in);
      if(name2proposal.find(name) == name2proposal.end())
	{
	  std::cerr << "Could not parse the checkpoint file.  A reason for this may be that\n"
		    << "you used a different configuration or alignment file in combination\n"
		    << "with this checkpoint file. Fatality." << std::endl; 
	  ParallelSetup::genericExit(-1); 
	}
      name2proposal[name]->readFromCheckpoint(in);
      ++ctr; 
    }  
}


void Chain::readFromCheckpoint( std::istream &in ) 
{
  chainRand.readFromCheckpoint(in); 
  couplingId = cRead<int>(in);   
  likelihood = cRead<double>(in); 
  lnPr = cRead<double>(in);   
  currentGeneration = cRead<int>(in); 

  initProposalsFromStream(in);

  std::unordered_map<std::string, AbstractParameter*> name2parameter; 
  for(auto &p : extractParameters())
    {
      std::stringstream ss;       
      p->printShort(ss);
      assert(name2parameter.find(ss.str()) == name2parameter.end()); // not yet there 
      name2parameter[ss.str()] = p;
    }  

  nat ctr = 0; 
  while(ctr < name2parameter.size())
    {
      std::string name = readString(in); 
      if(name2parameter.find(name) == name2parameter.end())
	{
	  std::cerr << "Could not parse the checkpoint file. A reason for this may be that\n"
		    << "you used a different configuration or alignment file in combination\n"
		    << "with this checkpoint file. Fatality." << std::endl; 
	  ParallelSetup::genericExit(-1); 
	}
      auto param  = name2parameter[name]; 

      ParameterContent content = param->extractParameter(*tralnPtr); // initializes the object correctly. the object must "know" how many values are to be extracted 
      content.readFromCheckpoint(in);
      savedContent[param->getId()]  = content; 

      ++ctr;
    }
}


void Chain::writeToCheckpoint( std::ostream &out) const
{
  chainRand.writeToCheckpoint(out);   
  cWrite(out, couplingId); 
  cWrite(out, likelihood); 
  cWrite(out, lnPr); 
  cWrite(out, currentGeneration); 

  for(auto &p : proposals)
    p->writeToCheckpoint(out);

  for(auto &var: extractParameters())
    {
      const auto &compo = savedContent.at(var->getId()); 

      std::stringstream ss; 
      var->printShort(ss); 

      std::string name = ss.str(); 

      writeString(out,name); 
      compo.writeToCheckpoint(out);
    }  
}   

// c++<3
