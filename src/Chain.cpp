#include <sstream>
#include <map> 
#include <unordered_map>

#include "proposals/StatNNI.hpp"
#include "proposals/NodeSlider.hpp"
#include "Chain.hpp"		
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"

#include "eval/LikelihoodEvaluator.hpp" 
#include "ParallelSetup.hpp"
#include "Category.hpp"
#include "TreePrinter.hpp"

// #define VERIFY_GEN 4
#define VERIFY_GEN 1000000


Chain:: Chain(randKey_t seed, std::shared_ptr<TreeAln> _traln, const std::vector<std::unique_ptr<AbstractProposal> > &_proposals, 
	      std::vector<ProposalSet> _proposalSets, LikelihoodEvaluator eval, bool isDryRun) 
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
  , evaluator(eval)
{
  for(auto &p : _proposals)
    {
      std::unique_ptr<AbstractProposal> copy(p->clone()); 
      assert(copy->getPrimaryParameterView().size() != 0); 
      proposals.push_back(std::move(copy)); 
    }

  auto root = tralnPtr->getAnyBranch();
  const std::vector<AbstractParameter*> vars = extractParameters(); 

  if(not isDryRun)
    evaluator.evaluate(*tralnPtr, root, true); 

  prior.initialize(*tralnPtr, vars);
  
  likelihood = tralnPtr->getTrHandle().likelihood; 
  suspend(); 
  updateProposalWeights();
}


Chain::Chain(  Chain&& rhs)   
  : tralnPtr(rhs.tralnPtr)
  , deltaT(rhs.deltaT)
  , runid(rhs.runid)
  , tuneFrequency(rhs.tuneFrequency)
  , hastings(rhs.hastings)
  , currentGeneration(rhs.currentGeneration)
  , couplingId(rhs.couplingId) 
  , proposals(std::move(rhs.proposals))
  , proposalSets(std::move(rhs.proposalSets))
  , chainRand(std::move(rhs.chainRand)) 
  , relWeightSumSingle(rhs.relWeightSumSingle)
  , relWeightSumSets(rhs.relWeightSumSets) 
  , prior(std::move(rhs.prior))
  , bestState(rhs.bestState)
  , evaluator(std::move(rhs.evaluator))
  , likelihood(rhs.likelihood)
  , lnPr(rhs.lnPr)
  , savedContent(std::move(rhs.savedContent))
{
  // tout << "move constructing" << std::endl; 
}


Chain::Chain(  const Chain& rhs)   
  : tralnPtr(rhs.tralnPtr)
  , deltaT(rhs.deltaT)
  , runid(rhs.runid)
  , tuneFrequency(rhs.tuneFrequency)
  , hastings(rhs.hastings)
  , currentGeneration(rhs.currentGeneration)
  , couplingId(rhs.couplingId) 
  , proposalSets(std::move(rhs.proposalSets))
  , chainRand(std::move(rhs.chainRand)) 
  , relWeightSumSingle(rhs.relWeightSumSingle)
  , relWeightSumSets(rhs.relWeightSumSets) 
  , prior(std::move(rhs.prior))
  , bestState(rhs.bestState)
  , evaluator(std::move(rhs.evaluator))
  , likelihood(rhs.likelihood)
  , lnPr(rhs.lnPr)
  , savedContent(std::move(rhs.savedContent))
{
  for(auto &p : rhs.proposals)
    proposals.emplace_back(p->clone());

  // tout << "copy constructing" << std::endl; 
}




Chain& Chain::operator=(Chain rhs)
{
  swap(*this, rhs); 
  return *this; 
}


void swap(Chain &lhs, Chain &rhs)
{
  using std::swap; 

  assert(lhs.currentGeneration == rhs.currentGeneration); 

  swap(lhs.tralnPtr,rhs.tralnPtr); 
  swap(lhs.deltaT,rhs.deltaT);
  swap(lhs.runid,rhs.runid);
  swap(lhs.tuneFrequency,rhs.tuneFrequency);
  swap(lhs.hastings,rhs.hastings);
  swap(lhs.currentGeneration,rhs.currentGeneration);
  swap(lhs.couplingId,rhs.couplingId); 
  swap(lhs.proposals,rhs.proposals);
  swap(lhs.proposalSets,rhs.proposalSets);
  swap(lhs.chainRand,rhs.chainRand);
  swap(lhs.relWeightSumSingle,rhs.relWeightSumSingle);
  swap(lhs.relWeightSumSets,rhs.relWeightSumSets); 
  swap(lhs.prior,rhs.prior);
  swap(lhs.bestState,rhs.bestState);
  swap(lhs.evaluator,rhs.evaluator);
  swap(lhs.likelihood,rhs.likelihood);
  swap(lhs.lnPr,rhs.lnPr);
  swap(lhs.savedContent,rhs.savedContent);
}


void swapHeatAndProposals(Chain &lhs, Chain& rhs)
{
  std::swap(lhs.couplingId, rhs.couplingId); 
  std::swap(lhs.proposalSets, rhs.proposalSets); 
  std::swap(lhs.proposals, rhs.proposals); 
}

std::ostream& Chain::addChainInfo(std::ostream &out) const 
{
  return out << "[run=" << runid << ",heat=" << couplingId << ",gen=" << getGeneration() << "]" ; 
}

 
double Chain::getChainHeat() const 
{
  double tmp = 1. + ( deltaT * couplingId ) ; 
  double inverseHeat = 1.f / tmp; 
  assert(couplingId == 0 || inverseHeat < 1.); 
  return inverseHeat; 
}

void Chain::updateProposalWeights()
{
  // prepare proposals 
  relWeightSumSingle = 0; 
  for(auto& elem : proposals)
    relWeightSumSingle +=  elem->getRelativeWeight();
  
  relWeightSumSets = 0; 
  for(auto &set : proposalSets)
    relWeightSumSets += set.getRelativeWeight();
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
	  v->applyParameter(*tralnPtr, savedContent[v->getId()]); 
	}
    }

  // now deal with all the other parameters 
  for(auto &v : vs)
    {
      if(v->getCategory() != Category::TOPOLOGY)
	v->applyParameter(*tralnPtr, savedContent[v->getId()]); 
    }

  if(evaluate)
    {
      auto root = tralnPtr->getAnyBranch(); 
      evaluator.evaluate(*tralnPtr, root, true);

      if(fabs(likelihood - tralnPtr->getTrHandle().likelihood) >  ACCEPTED_LIKELIHOOD_EPS ) 
	{
	  addChainInfo(std::cerr); 
	  std::cerr << "While trying to resume chain: previous chain liklihood"
		    <<  " could not be exactly reproduced. Please report this issue." << std::endl; 
	  std::cerr << MAX_SCI_PRECISION << 
	    "prev=" << likelihood << "\tnow=" << tralnPtr->getTrHandle().likelihood << std::endl; 	  
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

  updateProposalWeights();
}


ProposalSet& Chain::drawProposalSet(Randomness &rand)
{
  double r = relWeightSumSets * rand.drawRandDouble01(); 

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


AbstractProposal& Chain::drawProposalFunction(Randomness &rand)
{ 
  double r = relWeightSumSingle * rand.drawRandDouble01();   

  for(auto& c : proposals)
    {
      double w = c->getRelativeWeight(); 
      if(r < w )
	{
	  // std::cout << "drawn " << c.get() << std::endl; 
	  return *c;
	}
      else 
	r -= w; 
    }

  assert(0); 
  return *(proposals[0]); 
}


std::string Chain::serializeConditionally( CommFlag commFlags)  const 
{
  auto &&ss = std::stringstream{}; 
  
  if(commFlags & CommFlag::PrintStat)
    {
      ss.write(reinterpret_cast<const char*>(&couplingId), sizeof(couplingId)); 
      ss.write(reinterpret_cast<const char*>(&bestState), sizeof(bestState )); 
      ss.write(reinterpret_cast<const char*>(&likelihood), sizeof(likelihood)); 
      ss.write(reinterpret_cast<const char*>(&lnPr), sizeof(lnPr)); 
      ss.write(reinterpret_cast<const char*>(&currentGeneration), sizeof(currentGeneration)); 
    }

  if(commFlags & CommFlag::Proposals)
    {
      for(auto& p : proposals)
	p->serialize(ss); 

      for(auto& p: proposalSets)
	p.serialize(ss);
    }

  if(commFlags & CommFlag::Tree)
    {
      for(auto &var: extractParameters())
	{
	  const ParameterContent& compo = savedContent.at(var->getId()); 
	  
	  std::stringstream tmp; 
	  var->printShort(tmp); 	  	  
	  this->writeString(ss,tmp.str()); 
	  compo.serialize(ss);
	}
    }

  return ss.str(); 
}


void Chain::deserializeConditionally(std::string str, CommFlag commFlags)
{    
  auto &&ss = std::stringstream{}; 
  ss.str(str); 

  if(commFlags & CommFlag::PrintStat)
    {
      ss.read(reinterpret_cast<char*>(&couplingId), sizeof(couplingId) ); 
      ss.read(reinterpret_cast<char*>(&bestState), sizeof(bestState)); 
      ss.read(reinterpret_cast<char*>(&likelihood), sizeof(likelihood)); 
      ss.read(reinterpret_cast<char*>(&lnPr), sizeof(lnPr)); 
      ss.read(reinterpret_cast<char*>(&currentGeneration), sizeof(currentGeneration)); 
    }

  if(commFlags & CommFlag::Proposals)
    {
      initProposalsFromStream(ss);
    }
  
  if(commFlags & CommFlag::Tree)
    {
      auto name2parameter = std::unordered_map<std::string, AbstractParameter*>{} ; 
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
	  content.deserialize(ss);
	  savedContent[param->getId()]  = content; 

	  ++ctr;
	}
    }
}




BranchPlain Chain::peekNextVirtualRoot(TreeAln &traln, Randomness rand)  
{
  nat curGen = rand.getGeneration();
  // tout << rand << std::endl; 

  rand.rebaseForGeneration(curGen + 1 ); 

  // tout << rand << std::endl; 

  double sum = relWeightSumSingle + relWeightSumSets; 
  
  auto branch = BranchPlain(); 

  // return traln.getAnyBranch();

  if( rand.drawRandDouble01() * sum <  relWeightSumSingle) 
    {
      auto &pfun = drawProposalFunction(rand); 
      branch = pfun.determinePrimeBranch(traln, rand); 
    }
  else 
    {
      auto pset = drawProposalSet(rand); 
      branch = pset.getProposalView()[0]->determinePrimeBranch(traln, rand); 
    }

  if(branch.getPrimNode() == 0 || branch.getSecNode() == 0)
    {
      // TODO better draw  
      branch = TreeRandomizer::drawInnerBranchUniform(traln, rand); 
#ifdef PRINT_EVAL_CHOICE
      tout << "RANDOM " << branch<< std::endl; 
#endif
      // tout << "no prediction, just returning a somewhat inner branch " << branch  << std::endl; 
    }
  else 
    {
#ifdef PRINT_EVAL_CHOICE
      tout << "PREDICT " << branch << std::endl; 
#endif
      // tout << "predicted " << branch << std::endl;       
    }
  

  return branch; 
}



void Chain::stepSingleProposal()
{
  auto &traln = *tralnPtr; 
  
  double prevLnl = likelihood ; 
  double myHeat = getChainHeat();

  auto& pfun = drawProposalFunction(chainRand);

  /* reset proposal ratio  */
  hastings = 0; 

  pfun.applyToState(*tralnPtr, prior, hastings, chainRand, evaluator);


  // tout << " have "  << pfun << std::endl; 

  // tout << chainRand << std::endl; 
  auto suggestion = peekNextVirtualRoot(traln,chainRand); 
  // auto suggestion = traln.getAnyBranch() ; 
  // assert(0);  

  pfun.evaluateProposal(evaluator, *tralnPtr, suggestion);
  
  double priorRatio = prior.getLnPriorRatio();
  double lnlRatio = traln.getTrHandle().likelihood - prevLnl; 

  double testr = chainRand.drawRandDouble01();
  double acceptance = exp(( priorRatio + lnlRatio) * myHeat + hastings) ; 

  bool wasAccepted  = testr < acceptance; 

#ifdef DEBUG_SHOW_EACH_PROPOSAL 
  auto& output = tout  ; 
// std::cout ; 

  addChainInfo(output); 
  output << "\t" << (wasAccepted ? "ACC" : "rej" )  << "\t"<< pfun.getName() << "\t" 
	 << MORE_FIXED_PRECISION << prevLnl << "\tdelta(lnl)=" << lnlRatio << "\tdelta(lnPr)=" << priorRatio << "\thastings=" << hastings << std::endl; 
#endif

  if(wasAccepted)
    {
      pfun.accept();      
      prior.accept();
      if(bestState < traln.getTrHandle().likelihood  )
	bestState = traln.getTrHandle().likelihood; 
      likelihood = traln.getTrHandle().likelihood; 
      lnPr = prior.getLnPrior();
    }
  else
    {
      pfun.resetState(traln);
      pfun.reject();
      prior.reject();
      
      auto myRejected = std::vector<bool>(traln.getNumberOfPartitions(), false); 
      for(auto &elem : pfun.getAffectedPartitions())
	myRejected[elem] = true; 

      auto nodes = pfun.getInvalidatedNodes(traln); 
      evaluator.accountForRejection(traln, myRejected, nodes); 
    }
  
  // tout << "ORIENT at end: " << evaluator.getOrientation() << std::endl; 
  evaluator.freeMemory();


  if(this->tuneFrequency <  pfun.getNumCallSinceTuning() ) 
    pfun.autotune();
}


void Chain::stepSetProposal()
{
  auto& traln = *tralnPtr; 

  double myHeat = getChainHeat(); 
  auto &pSet = drawProposalSet(chainRand); 

  // tout << "have " << pSet << std::endl; 

  auto oldPartitionLnls = traln.getPartitionLnls(); 
  
  auto p2Hastings =  std::unordered_map<AbstractProposal*, double>{} ; 
  auto p2LnPriorRatio = std::unordered_map<AbstractProposal*, double>{}; 
  auto p2OldLnl = std::unordered_map<AbstractProposal*, double>{} ; 

  auto affectedPartitions = std::vector<nat>{}; 
  
  auto branches = pSet.getProposalView()[0]->prepareForSetExecution(traln, chainRand);

  for(auto &proposal : pSet.getProposalView())
    {
      proposal->setPreparedBranch(branches.first);
      proposal->setOtherPreparedBranch(branches.second);

      double lHast = 0;       
      prior.reject();
      proposal->applyToState(traln, prior, lHast, chainRand, evaluator); 
      p2LnPriorRatio[proposal] = prior.getLnPriorRatio(); 
      p2Hastings[proposal] = lHast; 

      for(auto &p : proposal->getPrimaryParameterView())
	{
	  auto partitions = p->getPartitions(); 
	  affectedPartitions.insert(end(affectedPartitions), begin(partitions), end(partitions)); 
	  
	  double lnl = 0; 
	  for(auto &partition : partitions)
	    lnl += oldPartitionLnls[partition]; 
	  p2OldLnl[proposal] = lnl; 
	}
    }

  bool fullTraversalNecessary = pSet.needsFullTraversal();

  if( branches.first.equalsUndirected(BranchPlain(0,0)) ) // TODO another HACK
    {
      // this should be a reasonable suggestion 
      auto nextRoot = peekNextVirtualRoot(traln, chainRand); 
      evaluator.evaluatePartitionsWithRoot(traln, nextRoot, affectedPartitions, fullTraversalNecessary); 
    }
  else 
    {
      pSet.getProposalView()[0]->prepareForSetEvaluation(traln, evaluator);
      evaluator.evaluatePartitionsWithRoot(traln, branches.first , affectedPartitions, fullTraversalNecessary );
    }

  auto newPLnls = traln.getPartitionLnls();

  prior.reject();		// slight abuse 
  auto partitionsToReset = std::vector<bool>(traln.getNumberOfPartitions() , false); 
  nat accCtr = 0; 
  nat total = 0; 
  auto p2WasAccepted = std::unordered_map<AbstractProposal*, bool>{} ; 
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
	  p2WasAccepted[proposal] = true; 
	}
      else 
	{
	  // TODO prior more efficient 
	  proposal->resetState(traln);
	  proposal->reject();	  

	  for(auto& param: proposal->getPrimaryParameterView())
	    {
	      for(auto p : param->getPartitions())
		partitionsToReset[p] = true; 
	    }
	  p2WasAccepted[proposal] = false; 
	}      
    }

  auto nodes = pSet.getProposalView()[0]->getInvalidatedNodes(traln); 
  evaluator.accountForRejection(traln, partitionsToReset, nodes); 

  likelihood = traln.getTrHandle().likelihood; 
  lnPr = prior.getLnPrior();
  
  if(bestState < likelihood)
    bestState = likelihood; 

#ifdef DEBUG_SHOW_EACH_PROPOSAL
  auto &output = std::cout ; 

  addChainInfo(output);
  output << "\t" << accCtr << "/"  << total << "\t" << pSet << "\t" << likelihood << std::endl; 
#endif

  for(auto &proposal : pSet.getProposalView())
    if(this->tuneFrequency < proposal->getNumCallSinceTuning()) // meh
      proposal->autotune(); 

  prior.reject();
  tout << MAX_SCI_PRECISION; 
  for(auto &proposal : pSet.getProposalView())
    {
      if(p2WasAccepted[proposal])
	{
	  double tmp = p2LnPriorRatio.at(proposal); 
	  prior.addToRatio(tmp); 
	}
    }

  evaluator.freeMemory();
  prior.accept();
}


void Chain::step()
{
  auto &traln = *tralnPtr; 
  currentGeneration++; 

  if(INTEGRATION_GENERATION < currentGeneration )
    {
      // startIntegration = true; 
      // assert(0); 
    }

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(traln, extractParameters());
#endif

  evaluator.imprint(traln);
  // inform the rng that we produce random numbers for generation x  
  chainRand.rebaseForGeneration(currentGeneration);

  double sum = relWeightSumSingle + relWeightSumSets; 
  if(chainRand.drawRandDouble01() * sum < relWeightSumSingle)
    stepSingleProposal();
  else 
    stepSetProposal();

#ifdef DEBUG_LNL_VERIFY
  evaluator.expensiveVerify(traln, traln.getAnyBranch() , likelihood); 
#endif

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(traln, extractParameters());
#endif

  // tout << "================================================================" << std::endl; 
}


void Chain::suspend()  
{
  auto params = extractParameters();
  savedContent.clear(); 

  for(auto& v : params)
    {
      assert(savedContent.find(v->getId()) == savedContent.end()); 
      savedContent[v->getId()] = v->extractParameter(*tralnPtr); 
    }
  
  tralnPtr->clearMemory(); 
  lnPr = prior.getLnPrior();
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
  auto result = std::unordered_set<AbstractParameter*>{}; 
  
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

  auto result2 = std::vector<AbstractParameter*> {}; 
  for(auto &v : result) 
    result2.push_back(v); 

  return result2; 
}


const std::vector<AbstractProposal*> Chain::getProposalView() const 
{
  auto result = std::vector<AbstractProposal*>{};  
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


void Chain::sample(  std::unordered_map<nat,TopologyFile> &paramId2TopFile ,  ParameterFile &pFile  )  const
{
  std::vector<AbstractParameter*> blParams; 
  for(auto &param : extractParameters()) 
    if(param->getCategory() == Category::BRANCH_LENGTHS)
      blParams.push_back(param); 

  for(auto &param : blParams)
    {
      nat myId = param->getIdOfMyKind(); 
      TopologyFile &f = paramId2TopFile.at(myId); 
      f.sample(*tralnPtr, getGeneration(), param); 
    }

  pFile.sample( *tralnPtr, extractParameters(), getGeneration(), prior.getLnPrior()); 
}



void Chain::initProposalsFromStream(std::istream& in)
{
  for(auto &p : proposals)
    {
      nat elem = cRead<int>(in); 
      assert(p->getId() == elem); 
      p->deserialize(in);
    }

  for(auto &p : proposalSets)
    p.deserialize(in);

  // auto name2proposal = std::unordered_map<std::string, AbstractProposal*>{}; 
  // for(auto &p :proposals)
  //   {
  //     std::stringstream ss; 
  //     p->printShort(ss); 
  //     assert(name2proposal.find(ss.str()) == name2proposal.end()); // not yet there 
  //     name2proposal[ss.str()] = p.get();
  //   }

  // nat ctr = 0; 
  // while(ctr < proposals.size())
  //   {
  //     std::string name = readString(in);
  //     if(name2proposal.find(name) == name2proposal.end())
  // 	{
  // 	  std::cerr << "Could not parse the checkpoint file.  A reason for this may be that\n"
  // 		    << "you used a different configuration or alignment file in combination\n"
  // 		    << "with this checkpoint file. Fatality." << std::endl; 
  // 	  ParallelSetup::genericExit(-1); 
  // 	}
  //     name2proposal[name]->readFromCheckpoint(in);
  //     ++ctr; 
  //   }
  
}


void Chain::deserialize( std::istream &in ) 
{
#ifdef DEBUG_SERIALIZE
  tout << "deserializing chain" << std::endl; 
#endif

  chainRand.deserialize(in); 
  couplingId = cRead<int>(in);   
  likelihood = cRead<double>(in); 
  lnPr = cRead<double>(in);   
  currentGeneration = cRead<int>(in); 

  initProposalsFromStream(in);

  auto name2parameter = std::unordered_map<std::string, AbstractParameter*> {}; 
  for(auto &p : extractParameters())
    {
      auto &&ss = std::stringstream{} ; 
      p->printShort(ss);
      assert(name2parameter.find(ss.str()) == name2parameter.end()); // not yet there 
      name2parameter[ss.str()] = p;
    }

  nat ctr = 0; 
  while(ctr < name2parameter.size())
    {
      auto name = readString(in); 
      if(name2parameter.find(name) == name2parameter.end())
	{
	  std::cerr << "Could not parse the checkpoint file. A reason for this may be that\n"
		    << "you used a different configuration or alignment file in combination\n"
		    << "with this checkpoint file. Fatality." << std::endl; 
	  ParallelSetup::genericExit(-1); 
	}
      auto param  = name2parameter[name]; 
      
      auto content = param->extractParameter(*tralnPtr); // initializes the object correctly. the object must "know" how many values are to be extracted 
      content.deserialize(in);
      savedContent[param->getId()]  = content; 

      ++ctr;
    }
}


void Chain::serialize( std::ostream &out) const
{
  chainRand.serialize(out);   
  cWrite(out, couplingId); 
  cWrite(out, likelihood); 
  cWrite(out, lnPr); 
  cWrite(out, currentGeneration); 

  for(auto &p : proposals)
    p->serialize(out);

  for(auto &var: extractParameters())
    {
      const auto &compo = savedContent.at(var->getId()); 

      std::stringstream ss; 
      var->printShort(ss); 

      std::string name = ss.str(); 

      writeString(out,name); 
      compo.serialize(out);
    }  
}   

// c++<3

