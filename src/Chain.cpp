#include <sstream>

#include "Chain.hpp"
#include "Topology.hpp"
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
  , state(*traln)
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
  , state(rhs.state)
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
  assert(traln->getNumBranches() == 1); 
  state.accessTopology().restoreTopology(*traln); 

  /* restore model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      Partition& partInfo = state.accessPartition(i);
      double alpha = partInfo.getAlpha(); 
      traln->setAlphaBounded(alpha,i) ; 
      
      std::vector<double> revMat = state.accessPartition(i).getRevMat();       
      traln->setRevMatBounded(revMat, i); 
      
      std::vector<double> stateFreqs = state.accessPartition(i).getStateFreqs(); 
      traln->setFrequenciesBounded(stateFreqs, i); 

      traln->initRevMat(i);
      traln->discretizeGamma(i);	 
    } 
  
  evaluator->evaluateFullNoBackup(*traln); 
  
  if(fabs(likelihood - traln->getTr()->likelihood) > ACCEPTED_LIKELIHOOD_EPS)
    {
      std::cerr << "While trying to resume chain: previous chain liklihood larger than " <<
	"evaluated likelihood. This is a programming error." << std::endl; 
      assert(0);       
    }  

  auto vs = extractVariables(); 
  prior.initialize(*traln, vs);

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
       << setprecision(2) << fixed << lnl  << endl; //   
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



void Chain::printParams(FILE *fh)
{
  tree *tr = traln->getTr(); 

  double treeLength = Branch(0,0,traln->getTreeLengthExpensive()).getInterpretedLength(*traln); 

  assert(treeLength != 0.); 
  fprintf(fh, "%d\t%f\t%f\t%.3f", currentGeneration,
	  prior.getLnPrior(), 
	  tr->likelihood,  
	  treeLength); 

  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      fprintf(fh, "\t%f\t%f", traln->accessPartitionLH( i),traln->getAlpha(i)) ; 

      auto revMat = traln->getRevMat(i);
      for(auto r : revMat)
	fprintf(fh, "\t%.2f", r); 
      
      auto freqs = traln->getFrequencies(i); 
      for(auto f : freqs)
	fprintf(fh, "\t%.2f", f); 
    }

  fprintf(fh, "\n"); 
  fflush(fh); 
}


void Chain::printParamFileStart(FILE *fh)
{
  fprintf(fh, "[ID: %s]\n", "TODO");
  fprintf(fh, "Gen\tLnPr\tLnL\tTL"); 

  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      fprintf(fh, "\tlnl.%d\talhpa.%d", i,i); 
      fprintf(fh, "\tr.%d(A<->C)\tr.%d(A<->G)\tr.%d(A<->T)\tr.%d(C<->G)\tr.%d(C<->T)\tr.%d(G<->T)",i, i,i,i,i,i); 
      fprintf(fh, "\tpi(A).%d\tpi(C).%d\tpi(G).%d\tpi(T).%d", i,i,i,i); 
    }

  fprintf(fh, "\n"); 
  fflush(fh); 
}


void Chain::finalizeOutputFiles(FILE *fh)
{
  /* topo file  */
  fprintf(fh, "end;\n"); 
}


void Chain::printSample(FILE *topofile, FILE *paramFile)
{
  this->printTopology(topofile);
  this->printParams(paramFile);
}


  /* TODO what about per model brach lengths? how does mrB do this? */
void Chain::printTopology(FILE *fh)
{  
  assert(couplingId == 0);
  
  TreePrinter tp(true, false, false);
  std::string treeString = tp.printTree(*traln);

  fprintf(fh,"\ttree gen.%d = [&U] %s\n", currentGeneration, treeString.c_str());
  fflush(fh);
}


void Chain::printNexusTreeFileStart( FILE *fh  )
{
  fprintf(fh, "#NEXUS\n[ID: %s]\n[Param: tree]\nbegin trees;\n\ttranslate\n", "TODO" ); 

  for(int i = 0; i < traln->getTr()->mxtips-1; ++i)
    fprintf(fh, "\t\t%d %s,\n", i+1, traln->getTr()->nameList[i+1]); 
  fprintf(fh, "\t\t%d %s;\n", traln->getTr()->mxtips, traln->getTr()->nameList[traln->getTr()->mxtips]);
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
  state.accessTopology().saveTopology(*(traln));

  /* save model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      Partition& partInfo =  state.accessPartition(i);       
      pInfo *partitionTr = traln->getPartition(i); 
      
      partInfo.setAlpha( partitionTr->alpha) ; 

      std::vector<double> tmp; 
      for(nat i = 0; i < numStateToNumInTriangleMatrix(partitionTr->states); ++i)
	tmp.push_back(partitionTr->substRates[i]); 
      partInfo.setRevMat(tmp); 
      tmp.clear(); 
      for(int i = 0; i < partitionTr->states ; ++i)
	tmp.push_back(partitionTr->frequencies[i]); 
      partInfo.setStateFreqs(tmp); 
    }  

  likelihood = traln->getTr()->likelihood; 
}

									   
std::vector<RandomVariable*> Chain::extractVariables() const 
{
  std::unordered_set<RandomVariable*> result; 
  for(auto &p : proposals)
    {
      for(auto &v : p->getPrimVar())
	result.insert(v); 
      for(auto &v : p->getSecVar())
	result.insert(v); 
    }
  
  std::vector<RandomVariable*> result2 ; 
  for(auto v : result)
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
  out << "\tLnL: " << rhs.getLnLikelihood() << "\tLnPr: " << rhs.getLnPrior(); 
  return out; 
}

