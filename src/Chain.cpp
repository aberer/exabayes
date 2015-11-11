#include <sstream>		
#include <map> 
#include <unordered_map>

// TODO remove 
#include "parameters/BranchLengthsParameter.hpp" 

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

#ifdef  _DEVEL
#include "memory.hpp" 
#endif

// #define VERIFY_GEN 4
#define VERIFY_GEN 1000000


Chain:: Chain(randKey_t seed, TreeAln traln, const std::vector<std::unique_ptr<AbstractProposal> > &proposals, 
	      std::vector<ProposalSet> proposalSets, LikelihoodEvaluator eval, bool isDryRun) 
  : _traln(traln)
  , _deltaT(0)
  , _runid(0)
  , _tuneFrequency(100)
  , _hastings(0)
  , _currentGeneration(0)
  , _couplingId(0)
  , _proposalSets(proposalSets)
  , _chainRand(seed)
  , _relWeightSumSingle(0)
  , _relWeightSumSets(0)
  , _bestState(std::numeric_limits<double>::lowest())
  , _evaluator(eval)
{
  for(auto &p : proposals)
    {
      std::unique_ptr<AbstractProposal> copy(p->clone()); 
      assert(copy->getPrimaryParameterView().size() != 0); 
      _proposals.push_back(std::move(copy)); 
    }

  const std::vector<AbstractParameter*> vars = extractParameters(); 
  _prior.initialize(_traln, vars);

  suspend(); 
  updateProposalWeights();
}


Chain::Chain(  Chain&& rhs)   
  : _traln(std::move(rhs._traln))
  , _deltaT(rhs._deltaT)
  , _runid(rhs._runid)
  , _tuneFrequency(rhs._tuneFrequency)
  , _hastings(rhs._hastings)
  , _currentGeneration(rhs._currentGeneration)
  , _couplingId(rhs._couplingId) 
  , _proposals(std::move(rhs._proposals))
  , _proposalSets(std::move(rhs._proposalSets))
  , _chainRand(std::move(rhs._chainRand)) 
  , _relWeightSumSingle(rhs._relWeightSumSingle)
  , _relWeightSumSets(rhs._relWeightSumSets) 
  , _prior(std::move(rhs._prior))
  , _bestState(rhs._bestState)
  , _evaluator(std::move(rhs._evaluator))
  , _lnPr(rhs._lnPr)
{
}


Chain::Chain(  const Chain& rhs)   
  : _traln(rhs._traln)
  , _deltaT(rhs._deltaT)
  , _runid(rhs._runid)
  , _tuneFrequency(rhs._tuneFrequency)
  , _hastings(rhs._hastings)
  , _currentGeneration(rhs._currentGeneration)
  , _couplingId(rhs._couplingId) 
  , _proposalSets(std::move(rhs._proposalSets))
  , _chainRand(std::move(rhs._chainRand)) 
  , _relWeightSumSingle(rhs._relWeightSumSingle)
  , _relWeightSumSets(rhs._relWeightSumSets) 
  , _prior(std::move(rhs._prior))
  , _bestState(rhs._bestState)
  , _evaluator(std::move(rhs._evaluator))
  , _lnPr(rhs._lnPr)
{
  for(auto &p : rhs._proposals)
    _proposals.emplace_back(p->clone());
}




Chain& Chain::operator=(Chain rhs)
{
  swap(*this, rhs); 
  return *this; 
}


void swap(Chain &lhs, Chain &rhs)
{
  using std::swap; 

  assert(lhs._currentGeneration == rhs._currentGeneration); 

  swap(lhs._traln,rhs._traln); 
  swap(lhs._deltaT,rhs._deltaT);
  swap(lhs._runid,rhs._runid);
  swap(lhs._tuneFrequency,rhs._tuneFrequency);
  swap(lhs._hastings,rhs._hastings);
  swap(lhs._currentGeneration,rhs._currentGeneration);
  swap(lhs._couplingId,rhs._couplingId); 
  swap(lhs._proposals,rhs._proposals);
  swap(lhs._proposalSets,rhs._proposalSets);
  swap(lhs._chainRand,rhs._chainRand);
  swap(lhs._relWeightSumSingle,rhs._relWeightSumSingle);
  swap(lhs._relWeightSumSets,rhs._relWeightSumSets); 
  swap(lhs._prior,rhs._prior);
  swap(lhs._bestState,rhs._bestState);
  swap(lhs._evaluator,rhs._evaluator);
  swap(lhs._lnPr,rhs._lnPr);
}


void swapHeatAndProposals(Chain &lhs, Chain& rhs)
{
  std::swap(lhs._bestState,    rhs._bestState); 
  std::swap(lhs._couplingId,   rhs._couplingId    ); 
  std::swap(lhs._proposals,    rhs._proposals     ); 
  std::swap(lhs._proposalSets, rhs._proposalSets  ); 
}


std::ostream& Chain::addChainInfo(std::ostream &out) const 
{
  return out << "[run=" << _runid << ",heat=" << _couplingId << ",gen=" << getGeneration() << "]" ; 
}

 
double Chain::getChainHeat() const 
{
  double tmp = 1. + ( _deltaT * _couplingId ) ; 
  double inverseHeat = 1.f / tmp; 
  assert(_couplingId == 0 || inverseHeat < 1.); 
  return inverseHeat; 
}


void Chain::updateProposalWeights()
{
  // prepare proposals 
  _relWeightSumSingle = 0; 
  for(auto& elem : _proposals)
    _relWeightSumSingle +=  elem->getRelativeWeight();
  
  _relWeightSumSets = 0; 
  for(auto &set : _proposalSets)
    _relWeightSumSets += set.getRelativeWeight();
}



void Chain::resume() 
{    
  auto vs = extractParameters(); 
  _prior.initialize(_traln, vs);
  double prNow = _prior.getLnPrior(); 
  
  if( fabs(_lnPr - prNow) > ACCEPTED_LNPR_EPS)
    {
      std::cerr << MAX_SCI_PRECISION ; 
      std::cerr << "While trying to resume chain: previous log prior could not be\n"
		<< "reproduced. This is a programming error." << std::endl; 
      std::cerr << "Prior was " << _lnPr << "\t now we have " << prNow << std::endl; 
      assert(0);       
    }

  updateProposalWeights();
}


ProposalSet& Chain::drawProposalSet(Randomness &rand)
{
  double r = _relWeightSumSets * rand.drawRandDouble01(); 

  for(auto &c : _proposalSets )
    {
      double w = c.getRelativeWeight(); 
      if(r < w )
	return c; 
      else 
	r -= w ; 
    }

  assert(0); 
  return _proposalSets[0]; 
}


AbstractProposal& Chain::drawProposalFunction(Randomness &rand)
{ 
  double r = _relWeightSumSingle * rand.drawRandDouble01();   

  for(auto& c : _proposals)
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
  return *(_proposals[0]); 
}


void Chain::serializeConditionally( std::ostream &out, CommFlag commFlags)  const 
{
  if( ( commFlags & CommFlag::PRINT_STAT )  != CommFlag::NOTHING)
    {
      cWrite<decltype(_couplingId)>(out, _couplingId);
      cWrite<decltype(_bestState)>(out, _bestState);

      double lnl = getLikelihood(); 
      cWrite<decltype(lnl)>(out, lnl); 

      cWrite<decltype(_lnPr)>(out, _lnPr); 
      cWrite<decltype(_currentGeneration)>(out,_currentGeneration);
    }

  if( (commFlags & CommFlag::RAND) != CommFlag::NOTHING)
    _chainRand.serialize(out);

  if( ( commFlags & CommFlag::PROPOSALS ) != CommFlag::NOTHING )
    {
      for(auto& p : _proposals)
	p->serialize(out); 

      for(auto& p: _proposalSets)
	p.serialize(out);
    }

  if( ( commFlags & CommFlag::TREE ) != CommFlag::NOTHING )
    {
      for(auto &var: extractSortedParameters())
	{
	  auto compo = var->extractParameter(_traln);
	  auto&&  tmp = std::stringstream{}; 
	  var->printShort(tmp); 	  	  
	  writeString(out,tmp.str()); 
	  compo.serialize(out);
	}
    }
}


void Chain::deserializeConditionally(std::istream& in, CommFlag commFlags)
{    
  if( ( commFlags & CommFlag::PRINT_STAT ) != CommFlag::NOTHING )
    {
      _couplingId = cRead<decltype(_couplingId)>(in); 
      _bestState = cRead<decltype(_bestState)>(in); 

      double lnl = cRead<decltype(lnl)>(in); 
      setLikelihood(lnl); 

      _lnPr = cRead<decltype(_lnPr)>(in); 
      _currentGeneration = cRead<decltype(_currentGeneration)>(in);
      
      // tout << "deser: " 
      // 	   << _couplingId << "\t"
      // 	   << _bestState << "\t"
      // 	   << lnl << "\t"
      // 	   << _lnPr << "\t"
      // 	   << _currentGeneration << std::endl ; 
    }

  if( (commFlags & CommFlag::RAND) != CommFlag::NOTHING)
    _chainRand.deserialize(in);

  if( ( commFlags & CommFlag::PROPOSALS ) != CommFlag::NOTHING )
    initProposalsFromStream(in);

  if( ( commFlags & CommFlag::TREE ) != CommFlag::NOTHING )
    {
      for(auto &param : extractSortedParameters())
	{
	  auto name = readString(in);
	  {
	    auto&& tmp = std::stringstream{}; 
	    param->printShort(tmp);
	    assert( tmp.str().compare(name) == 0 ); 
	  }
	  auto content = param->extractParameter(_traln); // initializes the object correctly. the object must "know" how many values are to be extracted 
	  content.deserialize(in);
	  param->applyParameter(_traln, content);
	}
    }
}




BranchPlain Chain::peekNextVirtualRoot(TreeAln &traln, Randomness rand)  
{
  // TODO go beyond next step  

  nat curGen = rand.getGeneration();
  // tout << rand << std::endl; 

  rand.rebaseForGeneration(curGen + 1 ); 

  double sum = _relWeightSumSingle + _relWeightSumSets; 
  
  auto branch = BranchPlain(); 

  if( rand.drawRandDouble01() * sum < _relWeightSumSingle) 
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
    }
  

  return branch; 
}



void Chain::stepSingleProposal()
{
  auto &traln = _traln; 
  
  // double prevLnl = _lnl ; 
  double prevLnl = getLikelihood(); 
  double myHeat = getChainHeat();

  auto& pfun = drawProposalFunction(_chainRand);

  // tout << "have "  << pfun << std::endl; 

  /* reset proposal ratio  */
  _hastings = 0; 

  pfun.applyToState(_traln, _prior, _hastings, _chainRand, _evaluator);

  auto suggestion = peekNextVirtualRoot(traln,_chainRand); 

  pfun.evaluateProposal(_evaluator, _traln, suggestion);


  double priorRatio = _prior.getLnPriorRatio();
  double lnlRatio = traln.getTrHandle().likelihood - prevLnl; 

  double testr = _chainRand.drawRandDouble01();
  double acceptance = exp(( priorRatio + lnlRatio) * myHeat + _hastings) ; 

  bool wasAccepted  = testr < acceptance; 

#ifdef DEBUG_SHOW_EACH_PROPOSAL 
  auto& output = tout  ; 
  addChainInfo(output); 
  output << "\t" << (wasAccepted ? "ACC" : "rej" )  << "\t"<< pfun.getName() << "\t" 
	 << MORE_FIXED_PRECISION << prevLnl << "\tdelta(lnl)=" << lnlRatio << "\tdelta(lnPr)=" << priorRatio << "\thastings=" << _hastings << std::endl; 
#endif

  if(wasAccepted)
    {
      pfun.accept();      
      _prior.accept();
      if(_bestState < traln.getTrHandle().likelihood  )
	_bestState = traln.getTrHandle().likelihood; 
      _lnPr = _prior.getLnPrior();
    }
  else
    {
      pfun.resetState(traln);
      pfun.reject();
      _prior.reject();
      
      auto myRejected = std::vector<bool>(traln.getNumberOfPartitions(), false); 
      for(auto &elem : pfun.getAffectedPartitions())
	myRejected[elem] = true; 

      auto nodes = pfun.getInvalidatedNodes(traln); 
      _evaluator.accountForRejection(traln, myRejected, nodes); 

    }


  _evaluator.freeMemory();

#ifdef _DEVEL
  auto totalMem = getCurrentMemory(); 

  nat used = 0, unused = 0; 
  std::tie(used, unused) =  _evaluator.getArrayReservoir().getUsedAndUnusedBytes(); 
  auto sum = used + unused; 

  auto overhead = totalMem - sum ; 

  tout << "MEM\t" << used / 1024  << "\t" << unused / 1024  << "\t" << overhead / 1024 << "\t" 
       << double(used) / double(totalMem)  << "\t"
       << double(unused) / double(totalMem)   << "\t"
       << double(overhead) / double(totalMem)  
       << "\t" << totalMem / 1024 
       << std::endl; 
#endif

  if(_tuneFrequency <  pfun.getNumCallSinceTuning() ) 
    pfun.autotune();
  
}


void Chain::stepSetProposal()
{
  auto& traln = _traln; 

  double myHeat = getChainHeat(); 
  auto &pSet = drawProposalSet(_chainRand); 

  // tout << "have " << pSet << std::endl; 

  auto oldPartitionLnls = traln.getPartitionLnls(); 
  
  auto p2Hastings =  std::unordered_map<AbstractProposal*, double>{} ; 
  auto p2LnPriorRatio = std::unordered_map<AbstractProposal*, double>{}; 
  auto p2OldLnl = std::unordered_map<AbstractProposal*, double>{} ; 

  auto affectedPartitions = std::vector<nat>{}; 
  
  auto branches = pSet.getProposalView()[0]->prepareForSetExecution(traln, _chainRand);

  for(auto &proposal : pSet.getProposalView())
    {
      proposal->setPreparedBranch(branches.first);
      proposal->setOtherPreparedBranch(branches.second);

      double lHast = 0;       
      _prior.reject();
      proposal->applyToState(traln, _prior, lHast, _chainRand, _evaluator); 
      p2LnPriorRatio[proposal] = _prior.getLnPriorRatio(); 
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
      auto nextRoot = peekNextVirtualRoot(traln, _chainRand); 
      _evaluator.evaluatePartitionsWithRoot(traln, nextRoot, affectedPartitions, fullTraversalNecessary, true); 
    }
  else 
    {
      pSet.getProposalView()[0]->prepareForSetEvaluation(traln, _evaluator);
      _evaluator.evaluatePartitionsWithRoot(traln, branches.first , affectedPartitions, fullTraversalNecessary, true );
    }

  auto newPLnls = traln.getPartitionLnls();

  _prior.reject();		// slight abuse 
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

      if(_chainRand.drawRandDouble01() < accRatio)
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
  _evaluator.accountForRejection(traln, partitionsToReset, nodes); 

  // _lnl = traln.getTrHandle().likelihood; 
  auto lnl = traln.getTrHandle().likelihood; 
  _lnPr = _prior.getLnPrior();
  
  if(_bestState < lnl)
    _bestState = lnl; 

#ifdef DEBUG_SHOW_EACH_PROPOSAL
  auto &output = std::cout ; 

  addChainInfo(output);
  output << "\t" << accCtr << "/"  << total << "\t" << pSet << "\t" << lnl << std::endl; 
#endif

  for(auto &proposal : pSet.getProposalView())
    if(this->_tuneFrequency < proposal->getNumCallSinceTuning()) // meh
      proposal->autotune(); 

  _prior.reject();
  tout << MAX_SCI_PRECISION; 
  for(auto &proposal : pSet.getProposalView())
    {
      if(p2WasAccepted[proposal])
	{
	  double tmp = p2LnPriorRatio.at(proposal); 
	  _prior.addToRatio(tmp); 
	}
    }

  _evaluator.freeMemory();
  _prior.accept();
}


void Chain::step()
{
  ++_currentGeneration; 

#ifdef DEBUG_VERIFY_LNPR
  _prior.verifyPrior(_traln, extractParameters());
#endif

  _evaluator.imprint(_traln);
  // inform the rng that we produce random numbers for generation x  
  _chainRand.rebaseForGeneration(_currentGeneration);

  double sum = _relWeightSumSingle + _relWeightSumSets; 
  if(_chainRand.drawRandDouble01() * sum < _relWeightSumSingle)
    stepSingleProposal();
  else 
    stepSetProposal();

#ifdef DEBUG_LNL_VERIFY
  _evaluator.expensiveVerify(_traln, _traln.getAnyBranch() , getLikelihood()); 
#endif

#ifdef DEBUG_VERIFY_LNPR
  _prior.verifyPrior(_traln, extractParameters());
#endif
}


void Chain::suspend()  
{
  _traln.clearMemory(_evaluator.getArrayReservoir()); 
  _lnPr = _prior.getLnPrior();
  _evaluator.freeMemory(); 
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


std::vector<AbstractParameter*> Chain::extractSortedParameters() const 
{
  auto result = extractParameters();
  std::sort(begin(result), end(result), 
	    [] (const AbstractParameter *lhs, const AbstractParameter* rhs)  -> bool 
	    {
	      auto prioLhs = lhs->getParamPriority();
	      auto prioRhs = rhs->getParamPriority();
	      if(prioLhs == prioRhs)
		return lhs->getId() < rhs->getId(); 
	      return prioLhs > prioRhs;  
	    });

  return result; 
}

std::vector<AbstractParameter*> Chain::extractParameters() const 
{
  auto result = std::unordered_set<AbstractParameter*>{}; 
  
  // add parameters in the default proposals 
  for(auto &p : _proposals)
    {
      for(auto &v : p->getPrimaryParameterView())
	result.insert(v); 
      for(auto &v : p->getSecondaryParameterView())
	result.insert(v); 
    }
  
  // add the proposals in the proposal sets 
  for(auto &p : _proposalSets)
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
  for(auto &elem: _proposals)
    result.push_back(elem.get()); 
  return result; 
}


std::ostream& operator<<(std::ostream& out, const Chain &rhs)
{
  rhs.addChainInfo(out); 
  // out  << "\tLnl: " << rhs.tralnPtr->getTr()->likelihood << "\tLnPr: " << rhs.prior.getLnPrior();
  out  << "\tLnl: " << rhs.getLikelihood() << "\tLnPr: " << rhs._lnPr; 
  return out;  
}


void Chain::sample(  std::unordered_map<nat,TopologyFile> &paramId2TopFile ,  ParameterFile &pFile ) const
{
  auto blParamsUnfixed=  std::vector<AbstractParameter*>{} ; 
  AbstractParameter *topoParamUnfixed = nullptr; 
  
  for(auto &param : extractParameters()) 
    {
      if(param->getCategory() == Category::BRANCH_LENGTHS && param->getPrior()->needsIntegration() )
	blParamsUnfixed.push_back(param); 
      else if(param->getCategory() == Category::TOPOLOGY && param->getPrior()->needsIntegration() )
	topoParamUnfixed = param  ; 
    }

  if(blParamsUnfixed.size() > 0)
    {
      for(auto &param : blParamsUnfixed)
	{
	  nat myId = param->getId(); 
	  auto &f = paramId2TopFile.at(myId); 
	  f.sample(_traln, getGeneration(), param); 
	}
    }
  else if(topoParamUnfixed != nullptr)
    {
      auto &f = paramId2TopFile.at(topoParamUnfixed->getId()); 
      f.sample(_traln, getGeneration(), topoParamUnfixed); 
    }

  pFile.sample( _traln, extractParameters(), getGeneration(), _prior.getLnPrior()); 
}



void Chain::initProposalsFromStream(std::istream& in)
{
  for(auto &p : _proposals)
    p->deserialize(in);

  for(auto &p : _proposalSets)
    p.deserialize(in);
}     


void Chain::serialize( std::ostream &out) const
{
  serializeConditionally(out, CommFlag::PRINT_STAT | CommFlag::PROPOSALS | CommFlag::TREE | CommFlag::RAND ) ;
}  


void Chain::deserialize( std::istream &in ) 
{
  deserializeConditionally(in, CommFlag::PRINT_STAT | CommFlag::PROPOSALS | CommFlag::TREE | CommFlag::RAND ) ;
}

// c++<3

