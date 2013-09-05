#include "StatNNI.hpp"
#include "Path.hpp"
#include "TreeRandomizer.hpp"
#include "Arithmetics.hpp"
#include "AdHocIntegrator.hpp"
#include "GibbsProposal.hpp"



// #define QUICK_HACK

#if defined(QUICK_HACK ) && not defined(_EXPERIMENTAL_INTEGRATION_MODE)
#error "need integration mode "
#endif

// #define _NNI_GIBBS
// #define _EXPERIMENTAL_GIBBS
// #define _EXPERIMENTAL_GIBBS_EXTRA

#if not defined(_EXPERIMENTAL_GIBBS) &&  defined(_EXPERIMENTAL_GIBBS_EXTRA)
#error "define _EXPERIMENTAL_GIBBS as well!"
#endif

// #define _EXPERIMENTAL_GIBBS

StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
#ifdef _EXPERIMENTAL_INTEGRATION_MODE
  , invocations(0)
#endif
{
  this->name = "stNNI" ; 
  this->category = Category::TOPOLOGY; 
  relativeWeight = 5; 
  needsFullTraversal = false; 
}

#ifdef _EXPERIMENTAL_INTEGRATION_MODE

#define LOCALMAP std::unordered_map<Branch, std::pair<double,double>, BranchHashNoLength, BranchEqualNoLength>  



static LOCALMAP
getStatisticsInEnvironment( const Branch &refBranch, const TreeAln &traln, nat numStep, nat thinning, nat depth, const std::vector<AbstractParameter*> params )
{
  std::unordered_map<Branch, std::pair<double,double>, BranchHashNoLength, BranchEqualNoLength> map; 

  std::vector<Branch> branches; 
  branches.push_back(refBranch); 
  // tout << "start: " << refBranch << std::endl; 
  std::vector<Branch> nextGuys = {refBranch, refBranch.getInverted()}; 
  for(nat i = 0; i < depth; ++i)
    {
      auto copy = nextGuys; 
      nextGuys.clear(); 
      for(auto &b : copy)
	{
	  if(not b.isTipBranch(traln))
	    {
	      auto desc = traln.getDescendents(b); 	  
	      nextGuys.push_back(desc.first.getInverted());
	      nextGuys.push_back(desc.second.getInverted());
	    }
	}
      branches.insert(branches.end(), nextGuys.begin(), nextGuys.end()); 
    }

  for(auto &b : branches)
    {
      b = traln.getBranch(b, params);
    }

  for(auto &b : branches)    
    {
      auto samples = ahInt->integrate(b,traln, numStep, thinning); 
      assert(map.find(b) == map.end()); 
      map[b] = Arithmetics::getMeanAndVar(samples); 
    }

// #endif

void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{    
  auto params = getBranchLengthsParameterView(); 
  assert(params.size( )== 1) ; 
  auto param = params[0]; 

  auto b = traln.getBranch(TreeRandomizer::drawInnerBranchUniform(traln, rand), param); 
  nodeptr p = b.findNodePtr(traln); 
  auto switchingBranch = BranchPlain( rand.drawRandDouble01() < 0.5  
				   ? p->back->next->back->number
				   : p->back->next->next->back->number, 
						      p->back->number ); 
  auto bls = std::vector<BranchPlain>{b.toPlain(), switchingBranch}; 
  move.extractMoveInfo(traln, bls, params); 

  std::vector<AbstractPrior* > priors; 
  bool multiplyBranches = false; 
  for(auto &v : secondaryParameters)
    {
      multiplyBranches |= v->getCategory() == Category::BRANCH_LENGTHS; 
      priors.push_back(v->getPrior()) ; 
#ifdef UNSURE
      // fixed prior? 
      assert(0); 
#endif
    }

#ifdef NO_SEC_BL_MULTI
  multiplyBranches = false; 
#endif

  if(multiplyBranches)
    {
      // move.multiplyBranches(traln, rand, hastings, prior, multiplier, priors); 
      assert(0);
    }

  // tout << sctr << std::endl; 

#ifdef QUICK_HACK
  auto branch = move.getEvalBranch(traln); 
  auto samples = ahInt->integrate(branch, traln, 10000, 10); 
  auto meanBefore =  Arithmetics::getMeanAndVar(samples).first; 
  auto lnlBefore = eval.evaluate(traln, traln.getAnyBranch(), true); 
  
  move.applyToTree(traln, params);
  
  samples  = ahInt->integrate(branch, traln, 10000, 10); 
  auto meanAfter = Arithmetics::getMeanAndVar(samples).first; 
  auto lnlAfter = eval.evaluate(traln, traln.getAnyBranch(), true); 

  tout << "NNI\t" << meanAfter / meanBefore << "\t" << lnlAfter - lnlBefore << std::endl; 
  
#else 
  move.applyToTree(traln, params);
#endif

#if defined( _EXPERIMENTAL_GIBBS) && defined(_NNI_GIBBS)
  nat numIter =20; 
  assert(params.size() == 1 ); 
  double afterMoveLnl = eval.evaluate(traln, b, true);

  auto newBranch = GibbsProposal::drawFromEsitmatedPosterior(b, eval, traln, rand, numIter, hastings,  params[0]); 
  traln.setBranch(newBranch, params[0]);

  prior.updateBranchLengthPrior(traln, b.getLength(params[0]), newBranch.getLength(params[0]), params[0]);
  
  // tout << "\t" << prevLnl << "\t" << afterMoveLnl << "\t" << eval.evaluate(traln,b,true) <<  "\t" << hastings << std::endl; 

#ifdef _EXPERIMENTAL_GIBBS_EXTRA
  auto desc = traln.getDescendents(b); 
  std::vector<Branch> bs ; 
  bs.push_back(desc.first); 
  bs.push_back(desc.second); 
  desc = traln.getDescendents(b.getInverted()); 
  bs.push_back(desc.first); 
  bs.push_back(desc.second); 
  
  for(auto &branch : bs)
    {
      branch = traln.getBranch(branch, params[0]); 
      eval.evaluate(traln, branch, false ); 
      auto newBranch = GibbsProposal::drawFromEsitmatedPosterior(branch, eval, traln, rand, numIter, hastings, params[0]); 
      traln.setBranch(newBranch, params[0]); 
      prior.updateBranchLengthPrior(traln , branch.getLength(params[0]), newBranch.getLength(params[0]), params[0]); 
    }

  // tout << "hastings " << hastings << std::endl; 
  LikelihoodEvaluator::disorientNode( b.findNodePtr(traln)); 
  LikelihoodEvaluator::disorientNode( b.getInverted().findNodePtr(traln)); 
#endif

#endif
}


void StatNNI::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln)
{
  auto evalBranch = move.getEvalBranch(traln); 
  nodeptr p = evalBranch.findNodePtr(traln);
  move.disorientAtNode(traln, p);
  evaluator.evaluate(traln, evalBranch, false); 
}

  evaluator.evaluate(traln, evalBranch, false); 
}

void StatNNI::resetState(TreeAln &traln)  
{
  move.revertTree(traln, getSecondaryParameterView()); 
  // debug_checkTreeConsistency(traln);
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
