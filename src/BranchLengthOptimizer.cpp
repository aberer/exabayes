#include "BranchLengthOptimizer.hpp"

#include "OptimizedParameter.hpp" 

#include "model/Branch.hpp" 
#include "comm/Communicator.hpp"

using std::vector; 


BranchLengthOptimizer::BranchLengthOptimizer(TreeAln& traln, const BranchPlain& branch, int maxIter, Communicator &comm, const std::vector<AbstractParameter*> &blParams )
  :  _comm(comm)
  , _branch(branch)
{
  for(auto i = 0u; i < blParams.size(); ++i)
    _optParams.emplace_back(traln, branch, blParams[i], maxIter);

  for(auto &param : blParams)
    _origBranches.push_back(traln.getBranch(branch,param));
}


bool BranchLengthOptimizer::hasConvergedAll()  const 
{
  auto result = true; 
  for(auto &v : _optParams )
    result &= v.hasFinished(); 
  return result; 
}


void BranchLengthOptimizer::applyToTraversalDescriptor(std::vector<bool> &execModel, TreeAln& traln) const 
{
  auto &tr = traln.getTrHandle(); 
  for(auto i = 0u ; i < execModel.size(); ++i )
    {
      auto& p = traln.getPartition(i); 
      tr.td[0].executeModel[i] = (execModel[i] && p.getWidth() > 0 )  ? PLL_TRUE : PLL_FALSE; 
    }

  for(auto &v : _optParams)
    v.applyValues(traln.getTrHandle().td[0].parameterValues);
}


void BranchLengthOptimizer::optimizeBranches(TreeAln &traln)  
{
  auto tr = &(traln.getTrHandle());
  auto pr = &(traln.getPartitionsHandle()); 

  int i,model; 
  bool firstIteration = true;

  auto dlnLdlz = std::vector<double>(traln.getNumberOfPartitions(),0 ); 
  auto d2lnLdlz2 = std::vector<double>(traln.getNumberOfPartitions(),0 ); 

  bool outerConverged = false; 
  while(not outerConverged)
    {
      auto execModel = std::vector<bool>(traln.getNumberOfPartitions(), false); 
      for(auto &v : _optParams)
	{
	  if(not v.hasFinished() )
	    {
	      if(v.isCurvatureOk())
		v.resetStep();
	      if(not v.isCurvatureOk())
		v.changeSide();
	    }
	}

      for(auto &v : _optParams)
	v.applyToMask(execModel);

#ifdef VERBOSE 
      tout << "executing: "; 
      for(auto v : execModel)
      	tout << (v ? 1 : 0) << ","; 
      tout << std::endl; 
#endif

      traln.setExecModel(execModel);
      applyToTraversalDescriptor(execModel, traln); 

      if(firstIteration)
	{
	  makenewzIterative(tr, pr);
	  firstIteration = false;
	}

      execCore(tr, pr, dlnLdlz.data(), d2lnLdlz2.data());

      dlnLdlz = _comm.get().allReduce(dlnLdlz);
      d2lnLdlz2 = _comm.get().allReduce(d2lnLdlz2); 

      for(auto &v : _optParams)
	{
	  if(not v.hasFinished())
	    {
	      v.extractDerivatives(traln, dlnLdlz, d2lnLdlz2);
	      if( not v.isCurvatureOk() )
		v.shortenBadBranch();
	    }
	}

      for(auto &oParam : _optParams)
	{
	  if(not oParam.hasFinished() && oParam.isCurvatureOk()  )
	    {
	      oParam.nrStep();
	      oParam.decrIter(); 
	      oParam.checkConvergence();
	    }
	}
      
      outerConverged = hasConvergedAll() ; 
    }
  
  // finally determine nrD1 and nrD2 at typical branch length for
  // those that are not converged
  // {
  //   auto execModel = std::vector<bool>(traln.getNumberOfPartitions(), false); 
  //   dlnLdlz = std::vector<double>(traln.getNumberOfPartitions(),0 );
  //   d2lnLdlz2 = std::vector<double>(traln.getNumberOfPartitions(),0 ); 
      
  //   for(auto &v : _optParams)
  //     {
  // 	// if(not v.hasConvergedNew())
  // 	if( v.getSecondDerivative() > 0 || v.getFirstDerivative() > 1 )
  // 	  {
  // 	    // tout << "NOT" << std::endl; 
  // 	    double typical = 0.1; 
	    
  // 	    v.setToTypicalBranch(typical, traln);
  // 	    v.applyToMask(execModel);
  // 	  }
  //     }

  //   applyToTraversalDescriptor(execModel, traln); 
  //   execCore(tr, pr, dlnLdlz.data(), d2lnLdlz2.data()); 

  //   for(auto &v : _optParams)
  //     {
  // 	if( v.getSecondDerivative() > 0 || v.getFirstDerivative() > 1 )
  // 	  {
  // 	    double prev = v.getFirstDerivative(); 
  // 	    v.extractDerivatives(traln, dlnLdlz, d2lnLdlz2); 
  // 	    tout << MAX_SCI_PRECISION << "len= 0.1\t" << v.getFirstDerivative() <<  "\t" << prev << "\t" << v.getSecondDerivative() << std::endl; 
  // 	  }
  //     }
  // }


  // reset stuff 
  traln.setExecModel( std::vector<bool>(traln.getNumberOfPartitions(), true ) ); 
  for(auto i = 0u; i < _blParams.size() ; ++i) 
    traln.setBranch(_origBranches[i], _blParams[i]); 
}



