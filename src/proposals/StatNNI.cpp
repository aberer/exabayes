#include "StatNNI.hpp"
#include "Path.hpp"
#include "TreeRandomizer.hpp"
#include "AdHocIntegrator.hpp"

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
      map[b] = ahInt->getMeanAndVar(samples); 
    }

  return map; 
}

static void printBeforeAndAfter(const Branch& refBranch, const TreeAln& traln, LOCALMAP &mapBefore , LOCALMAP &mapAfter, nat depth, const std::vector<AbstractParameter*> &params, nat invocations, bool accepted, std::ostream &out)
{ 
  assert(params.size() == 1); 

  auto before =  mapBefore.at(refBranch),
    after = mapAfter.at(refBranch); 

  out << "DEG-0\t"   <<  invocations << "\t" << (accepted ? "ACC" : "rej") << "\t"
       << traln.getBranch(refBranch, params[0]).getInterpretedLength(traln, params[0])
       << "\t" << before.first << "\t" << after.first
       << "\t" << before.second << "\t" << after.second << std::endl; 

  std::vector<Branch>current = { refBranch, refBranch.getInverted() };
  std::vector<Branch>next;
  for(auto &b : current)
    {
      auto desc = traln.getDescendents(b); 
      next.push_back(desc.first.getInverted()); 
      next.push_back(desc.second.getInverted()); 
    }
  current = next ; 
  next.clear(); 

  // first order 
  for(auto &b : current)
    {
      std::pair<double,double> 
	before, 
	after; 

      if(mapBefore.find(b) != mapBefore.end())
	before = mapBefore.at(b); 
      else
	{
	  auto tmp = b ; 
	  tmp.setSecNode(refBranch.getOtherNode(tmp.getSecNode()));
	  assert(mapBefore.find(tmp) != mapBefore.end());
	  before = mapBefore.at(tmp);
	}
      
      Branch afterBranch ;  

      if(mapAfter.find(b) != mapAfter.end())
	{
	  after = mapAfter.at(b);
	  afterBranch = b; 
	}
      else 
	{
	  auto tmp = b; 
	  tmp.setSecNode(refBranch.getOtherNode(tmp.getSecNode())); 
	  assert(mapAfter.find(tmp) != mapAfter.end()); 
	  after = mapAfter.at(tmp);
	  afterBranch = tmp; 
	}

      out << "DEG-1\t"   << invocations << "\t"  << (accepted ? "ACC" : "rej") << "\t"
	  << traln.getBranch(b,params[0]).getInterpretedLength(traln, params[0]) 
	  << "\t"  << before.first << "\t" << after.first
	  << "\t" << before.second << "\t" << after.second << std::endl; 

      if(not b.isTipBranch(traln))
	{
	  auto desc = traln.getDescendents(b);
	  // tout << desc.first << "\t" << desc.second << std::endl; 
	  next.push_back(desc.first.getInverted());
	  next.push_back(desc.second.getInverted());
	}
    }

  for(nat i = 2 ; i < depth+1; ++i)
    {
      current = next ; 
      next.clear();

      // second order  
      for(auto &b : current)
	{
	  auto before = mapBefore.at(b); 
	  auto after = mapAfter.at(b);
	 
	  out << "DEG-"  << i  << "\t"   << invocations << "\t" << (accepted ? "ACC" : "rej") << "\t"
	       << traln.getBranch(b,params[0]).getInterpretedLength(traln, params[0]) 
	       << "\t"  << before.first << "\t" << after.first
	       << "\t" << before.second << "\t" << after.second << std::endl;       

	  if(not b.isTipBranch(traln))
	    {
	      auto desc = traln.getDescendents(b); 
	      next.push_back(desc.first.getInverted());
	      next.push_back(desc.second.getInverted());
	    }
	}
    }
}



#endif

void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{    
  auto params = getBranchLengthsParameterView(); 

  Branch b =  TreeRandomizer::drawInnerBranchUniform(traln, rand) ; 
  b = traln.getBranch(b, params); 

  nodeptr p = b.findNodePtr(traln); 

  Branch switchingBranch = Branch( rand.drawRandDouble01() < 0.5  
				   ? p->back->next->back->number
				   : p->back->next->next->back->number, 
				   p->back->number ); 
  
  move.extractMoveInfo(traln, {b, switchingBranch}, params); 

#ifdef _EXPERIMENTAL_INTEGRATION_MODE
  nat numGen = 5e3; 
  nat step = 1e1; 
  nat depth = 3; 
  LOCALMAP samplesBefore; 
  if(startIntegration)
    samplesBefore = getStatisticsInEnvironment(b, traln, numGen, step, depth, params); 
  double prevLnl = traln.getTr()->likelihood; 
#endif

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

  
  move.applyToTree(traln, params);

#ifdef _EXPERIMENTAL_INTEGRATION_MODE
  // simulate acceptance rejection decision 

  bool accepted = ahInt->decideUponAcceptance(traln, prevLnl ); 

  if(startIntegration)
    {
      auto samplesAfter = getStatisticsInEnvironment(b, traln, numGen, step, depth, params); 
      tout << MAX_SCI_PRECISION; 
      printBeforeAndAfter(b, traln, samplesBefore, samplesAfter, depth, params, invocations, accepted, nniOut);
    }
  ++invocations; 
#endif


#ifdef _EXPERIMENTAL_GIBBS
  nat numIter =20; 
  assert(params.size() == 1 ); 
  eval.evaluate(traln, b, true);
  auto newBranch = GibbsProposal::drawFromEsitmatedPosterior(b, eval, traln, rand, numIter, hastings,  params[0]); 
  traln.setBranch(newBranch, params[0]);
  prior.updateBranchLengthPrior(traln, b.getLength(params[0]), newBranch.getLength(params[0]), params[0]);

  // auto desc = traln.getDescendents(b); 
  // std::vector<Branch> bs ; 
  // bs.push_back(desc.first); 
  // bs.push_back(desc.second); 
  // desc = traln.getDescendents(b.getInverted()); 
  // bs.push_back(desc.first); 
  // bs.push_back(desc.second); 
  
  // for(auto &branch : bs)
  //   {
  //     branch = traln.getBranch(branch, params[0]); 
  //     eval.evaluate(traln, branch, false ); 
  //     auto newBranch = GibbsProposal::drawFromEsitmatedPosterior(branch, eval, traln, rand, numIter, hastings, params[0]); 
  //     traln.setBranch(newBranch, params[0]); 
  //     prior.updateBranchLengthPrior(traln , branch.getLength(params[0]), newBranch.getLength(params[0]), params[0]); 
  //   }

  // tout << "hastings " << hastings << std::endl; 
  // LikelihoodEvaluator::disorientNode( b.findNodePtr(traln)); 
  // LikelihoodEvaluator::disorientNode( b.getInverted().findNodePtr(traln)); 
#endif
}


void StatNNI::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln)
{
  Branch evalBranch = move.getEvalBranch(traln); 
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
