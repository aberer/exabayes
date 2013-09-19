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


StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "stNNI" ; 
  this->category = Category::TOPOLOGY; 
  relativeWeight = 5; 
  needsFullTraversal = false; 
}


// #endif

void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{    
  auto params = getBranchLengthsParameterView(); 

  auto b = traln.getBranch(TreeRandomizer::drawInnerBranchUniform(traln, rand), params); 
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


void StatNNI::resetState(TreeAln &traln)  
{
  move.revertTree(traln, getSecondaryParameterView()); 
  // debug_checkTreeConsistency(traln);
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
