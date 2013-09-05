#include "TreeIntegrator.hpp"
#include "BranchLengthMultiplier.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "RestoringLnlEvaluator.hpp"
#include "Branch.hpp"
#include "priors/ExponentialPrior.hpp"
#include "ProposalRegistry.hpp"
#include "proposals/BranchLengthMultiplier.hpp"
#include "NniMove.hpp"
#include "Arithmetics.hpp"

TreeIntegrator::TreeIntegrator(std::shared_ptr<TreeAln> tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed)
{
  auto restorer = std::make_shared<ArrayRestorer>(*tralnPtr);  
  auto eval = std::unique_ptr<LikelihoodEvaluator>( new RestoringLnlEvaluator(restorer)); 

#ifdef DEBUG_LNL_VERIFY
  eval->setDebugTraln(debugTree);
#endif

  // s.t. we do not have to care about the branch length linking problem 
  assert(tralnPtr->getNumberOfPartitions() == 1 ); 


  std::vector<std::unique_ptr<AbstractParameter> > params; 
  params.emplace_back(std::unique_ptr<AbstractParameter>(new BranchLengthsParameter(0,0))); 
  for(nat i = 0; i < tralnPtr->getNumberOfPartitions(); ++i)
    params[0]->addPartition(i);



  double lambda = 10;
  params[0]->setPrior(std::make_shared< ExponentialPrior>(lambda));


  std::vector<std::unique_ptr<AbstractProposal> >  proposals;   
  proposals.emplace_back(new BranchLengthMultiplier(ProposalRegistry::initBranchLengthMultiplier)); 
  proposals[0]->addPrimaryParameter( std::move(params[0])); 

  std::vector<ProposalSet> pSets; 

  integrationChain = std::unique_ptr<Chain>( new Chain(seed, tralnPtr, proposals, pSets, std::move(eval) ));
  integrationChain->getEvaluator().imprint(*tralnPtr);  
}


// void TreeIntegrator::copyTree(const TreeAln &traln)
// {
//   TreeAln &myTree = integrationChain->getTraln();
//   myTree.copyModel(traln);
//   auto eval = integrationChain->getEvaluator();
//   auto tr = myTree.getTr(); 
//   eval->evaluate(myTree, Branch(tr->start->number,tr->start->back->number ), true );
// }

void TreeIntegrator::prepareChain( const TreeAln &otherTree, bool copyEverything)
{
  auto &traln = integrationChain->getTraln();
  if(copyEverything)
    traln.copyModel(otherTree); 
  
  integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(), true); 
  integrationChain->reinitPrior(); 
}




// tree must be prepared! 
Branch2Stat TreeIntegrator::integrateTree(nat numgen, nat thinning)
{
  auto samples = std::unordered_map<BranchLength, std::vector<double>>{}; 
  auto params = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  auto& traln = *(integrationChain->getTralnPtr()); 
  assert(params.size( )== 1); 
  auto param = params[0]; 

  for(nat i = 0; i < numgen; ++i)
    {
      integrationChain->step();

      if(i % thinning == 0)
	{
	  for(auto branch : integrationChain->getTralnPtr()->extractBranches(param))
	    {
	      branch = traln.getBranch(branch.toPlain(), param); 
	      samples[branch].push_back(branch.getInterpretedLength(traln,param)); 
	    }
	}
    }

  auto result = Branch2Stat{}; 
  for(auto sample : samples)  
    {
      result[sample.first] = make_pair(Arithmetics::getMean(sample.second), 
				       Arithmetics::getVariance(sample.second));
    }

  return result; 
}



#define NUM_GEN_BASE 10000


void TreeIntegrator::integrateAllBranches(const TreeAln &otherTree, std::string runid )
{
  auto&& ss = std::stringstream{}; 
  ss << PROGRAM_NAME << "_nniIntegration." << runid ; 
  std::ofstream outfile(ss.str());
  
  outfile << MAX_SCI_PRECISION; 

  prepareChain(otherTree,true); 
  auto& traln = *(integrationChain->getTralnPtr());
  auto params = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  assert(params.size() == 1 ) ;
  auto param = params[0]; 

  auto branch2id = std::unordered_map<BranchLength,nat>{};

  nat ctr = 0; 
  for(auto branch : traln.extractBranches(param))
    {
      branch2id[branch] = ctr; 
      ++ctr; 
    }

  auto numBr =  traln.getNumberOfBranches(); 
  auto numGen =  numBr * NUM_GEN_BASE;   
  auto thinning = numBr * 10 ; 

  prepareChain(otherTree, false); 
  auto initBranch2Stat = integrateTree(numGen , thinning);
  

  double initLnl = integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(),  true); 
  for(auto elem : initBranch2Stat)
    {
      auto branch = elem.first; 
      branch.setConvertedInternalLength( traln, param, elem.second.first );
      traln.setBranch(branch, param); 
    }
  double optLnl = integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(),  true); 
  

  
  // print 
  for(auto branch : traln.extractBranches(param))
    {
      outfile << branch2id.at(branch)
	      << "\t" << branch.toPlain(); 

      outfile << "\t0"  ; 
      
      auto stat = initBranch2Stat[branch]; 
      outfile << "\t" << stat.first ; 
      outfile << "\t" << stat.second ; 
      
      outfile << "\tINIT" << "\tINIT"
	      << "\t" << initLnl 
	      << "\t" << optLnl 
	      << std::endl; 
    }

  for(auto moveBranch : traln.extractBranches(param))
    {
      if( moveBranch.isTipBranch(traln))
	continue; 

      // the 2 NNI alternatives 
      auto backP = moveBranch.findNodePtr(traln)->back; 
      auto nniAlternativeBranches = { BranchPlain(backP->next->back->number, backP->number) ,
				      BranchPlain(backP->next->next->back->number, backP->number )} ; 

      bool isFirst = true; 
      for(auto nniAlternativeBranch : nniAlternativeBranches)
	{
	  auto move = NniMove{}; 
	  move.extractMoveInfo(traln, {moveBranch.toPlain(), nniAlternativeBranch.toPlain() } , params); 

	  move.applyToTree(traln, params);
	  prepareChain(otherTree, false); 

	  auto afterMove = integrateTree(numGen, thinning); 

	  double initLnl = integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(),  true); 
	  for(auto elem : afterMove)
	    {
	      auto branch = elem.first; 
	      branch.setConvertedInternalLength( traln, param, elem.second.first );
	      traln.setBranch(branch, param);
	    }
	  double optLnl = integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(),  true); 

	  for(auto b : traln.extractBranches(param))
	    {
	      auto b4 = move.mapToBranchBeforeMove(b); // that's not a bad variable name, but a star trek:nemesis  reference
	      auto dist = b.toPlain().getDistance(moveBranch.toPlain(), traln) ; 
	 
	      if(branch2id.find(b4.toBlDummy()) == branch2id.end())
		{
		  outfile << "could not find " << b4 << " mapped from " << b<< std::endl; 
		  assert(0); 
		}
     
	      outfile << branch2id.at(b4.toBlDummy()) 
		      << "\t" << b4.toPlain()
		      << "\t" << dist 
		      << "\t" << afterMove[b].first
		      << "\t" << afterMove[b].second
		      << "\t" << (isFirst ? "A" : "B") 
		      << "\t" << branch2id.at(moveBranch)
		      << "\t" << initLnl 
		      << "\t" << optLnl
		   << std::endl ; 

	    }
	  
	  move.revertTree(traln,params); 

	  isFirst = false; 
	}
    }
}

