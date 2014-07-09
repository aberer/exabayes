#include "DistributionBranchLength.hpp"
#include "system/BoundsChecker.hpp"
#include "priors/AbstractPrior.hpp"
#include "comm/ParallelSetup.hpp"
#include "DistributionProposer.hpp"
#include "BranchLengthOptimizer.hpp"

// ONLY gamma right now 


template<class C> 
DistributionBranchLength<C>::DistributionBranchLength( ) 
  : BranchLengthMultiplier(0)
  , _convTuner{1.61, 0.01,10, 0.1, false }
  , _nonConvTuner{   2. , 0.1, 10, 0.1 , false }	
{
  _name = std::string("blDist") + C::getName() ; 
  _category = Category::BRANCH_LENGTHS; 
  _usingOptimizedBranches = true; 
}
 


template<class C>
std::pair<BranchPlain,BranchPlain> DistributionBranchLength<C>::prepareForSetExecution(TreeAln &traln, Randomness &rand) 
{
  assert(_inSetExecution); 

  auto chosenBranch = determinePrimeBranch(traln,rand);
  auto dummyBranch = BranchPlain(0,0); 

  return std::make_pair(chosenBranch, dummyBranch);
} 


template<class C>
void DistributionBranchLength<C>::createProposer(TreeAln &traln, LikelihoodEvaluator& eval, BranchPlain b)
{
  auto params = getBranchLengthsParameterView();
  assert(params.size() == 1 );     
  auto param = params[0]; 
  
  _savedBranch = traln.getBranch(b.toPlain(),param); 
  
  eval.evaluateSubtrees(traln, b.toPlain(), param->getPartitions(), false);
  eval.evaluatePartitionsDry(traln, b.toPlain(), param->getPartitions());

  auto blo= BranchLengthOptimizer(traln, b.toPlain(), 30, eval.getParallelSetup().getChainComm(), params); 
  blo.optimizeBranches(traln);
  auto optParams = blo.getOptimizedParameters();
  assert(optParams.size() == 1);
  _proposer = optParams[0].getProposerDistribution< C >(traln, _convTuner.getParameter(), _nonConvTuner.getParameter());
}


template<class C>
void DistributionBranchLength<C>::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto b = proposeBranch(traln, rand).toBlDummy(); 
  
  auto params = getBranchLengthsParameterView();
  assert(params.size() == 1 );     
  auto param = params[0]; 

  _savedBranch = traln.getBranch(b.toPlain(), param);

  if(not _inSetExecution)
    createProposer(traln, eval, b.toPlain());

  auto newBranch  = _proposer.proposeBranch(b.toPlain(), traln, param, rand); 

  _converged = _proposer.isConverged(); 

  auto prevAbsLen = _savedBranch.getInterpretedLength(traln,param); 
  auto curAbsLen  = newBranch.getInterpretedLength(traln,param); 

  auto oldPr = param->getPrior()->getLogProb( ParameterContent{{ prevAbsLen }} ); 
  auto newPr = param->getPrior()->getLogProb( ParameterContent{{ curAbsLen }} ); 

  // update prior 
  prior.addToRatio( newPr / oldPr );

  // update hastings 
  auto propDensCur = _proposer.getLogProbability(curAbsLen); 
  auto propDensPrev = _proposer.getLogProbability(prevAbsLen);

  hastings *= (propDensPrev / propDensCur) ;

  traln.setBranch(newBranch, param);
}


template<class C>
void DistributionBranchLength<C>::autotune()
{
  auto res = getEnvironmentVariable("TUNE"); 

  if(res.compare("0") != 0 )
    {
      if( _convTuner.getRecentlySeen()  > 100 )
	{
	  _convTuner.tune();
	}
      
      if( _nonConvTuner.getRecentlySeen() > 100)
      	{
      	  double oldParam = _nonConvTuner.getParameter();
      	  double ratio = _nonConvTuner.getRatio(); 
      	  _nonConvTuner.tune();
      	  double newParam = _nonConvTuner.getParameter(); 
      	}
    }

  _sctr.nextBatch();
}

 

template<class C>
void DistributionBranchLength<C>::accept() 
{
  _sctr.accept();  

  if(_converged)
    _convTuner.accept();
  else 
    _nonConvTuner.accept();
}


template<class C>
void DistributionBranchLength<C>::reject() 
{
  _sctr.reject();

  if( _converged )
    _convTuner.reject();
  else 
    _nonConvTuner.reject();
}

