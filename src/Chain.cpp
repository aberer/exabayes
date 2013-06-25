#include <sstream>


#include "Chain.hpp"
#include "Topology.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "GlobalVariables.hpp"
#include "BipartitionHash.hpp"
#include "tune.h"
#include "ExtendedTBR.hpp"
#include "ExtendedSPR.hpp"
#include "ParsimonySPR.hpp" 
#include "TreeLengthMultiplier.hpp"
#include "StatNNI.hpp"
#include "PartitionProposal.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"
#include "LikelihoodEvaluator.hpp" 


// #define DEBUG_ACCEPTANCE

Chain::Chain(randKey_t seed, int id, int _runid, TreeAlnPtr _traln, const vector<ProposalPtr> &_proposals, int _tuneFreq ,const vector<RandomVariablePtr> &_variables, LikelihoodEvaluatorPtr eval) 
  : traln(_traln)
  , deltaT(0)
  , runid(_runid)
  , tuneFrequency(_tuneFreq)
  , hastings(0)
  , currentGeneration(0)
  , couplingId(id)
  , state(*traln)
  , chainRand(seed.v[0])
  , relWeightSum(0)
  , bestState(numeric_limits<double>::lowest())
  , prior(*_traln, _variables)
  , evaluator(eval)
  , variables(_variables)
{
  for(auto &p : _proposals)
    proposals.push_back(ProposalPtr(p->clone())); 

  for(int j = 0; j < traln->getNumberOfPartitions(); ++j)
    traln->initRevMat(j);

  saveTreeStateToChain(); 

  for(auto& elem : proposals)
    relWeightSum +=  elem->getRelativeWeight();

  cout << "relative weight is "  << relWeightSum << endl; 

  
  // we could do that with unique_ptr later ... 
  // for(auto &v:  _variables) 
  //   variables.push_back(v); 
}


ostream& Chain::addChainInfo(ostream &out)
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


void Chain::saveTreeStateToChain()
{
  state.accessTopology().saveTopology(*(traln));

  /* save model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      Partition& partInfo =  state.accessPartition(i);       
      pInfo *partitionTr = traln->getPartition(i); 
      
      partInfo.setAlpha( partitionTr->alpha) ; 

      vector<double> tmp; 
      for(int i = 0; i < numStateToNumInTriangleMatrix(partitionTr->states); ++i)
	tmp.push_back(partitionTr->substRates[i]); 
      partInfo.setRevMat(tmp); 
      tmp.clear(); 
      for(int i = 0; i < partitionTr->states ; ++i)
	tmp.push_back(partitionTr->frequencies[i]); 
      partInfo.setStateFreqs(tmp); 
    }  
}


void Chain::applyChainStateToTree()
{
  assert(traln->getNumBranches() == 1); 
  state.accessTopology().restoreTopology(*traln); 

  /* restore model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      Partition& partInfo = state.accessPartition(i);
      double alpha = partInfo.getAlpha(); 
      traln->setAlphaBounded(alpha,i) ; 
      
      vector<double> revMat = state.accessPartition(i).getRevMat();       
      traln->setRevMatBounded(revMat, i); 
      
      vector<double> stateFreqs = state.accessPartition(i).getStateFreqs(); 
      traln->setFrequenciesBounded(stateFreqs, i); 

      traln->initRevMat(i);
      traln->discretizeGamma(i);	 
    } 
  
  evaluator->evaluateFullNoBackup(*traln); 
  prior.reinitPrior(*traln, variables);
}




void Chain::debug_printAccRejc(ProposalPtr &prob, bool accepted, double lnl, double lnPr ) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
    tout << "[run=" << runid << ",heat="  << couplingId << ",gen="  << currentGeneration << "]\t" << (accepted ? "ACC" : "rej" )  << "\t"<< prob->getName() << "\t" << setprecision(2) << fixed << lnl  << endl; //   "\t"<< lnPr << endl; 
#endif
}




ProposalPtr& Chain::drawProposalFunction()
{ 
  double r = relWeightSum * chainRand.drawRandDouble01();   

  for(auto& c : proposals)
    {
      double w = c->getRelativeWeight(); 
      if(r < w )
	return c;
      else 
	r -= w; 
    }

  assert(0); 
  return proposals[0]; 
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
  string treeString = tp.printTree(*traln);

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
  swap(couplingId, rhs.couplingId); 
  swap(bestState, rhs.bestState); 
  // swap(chainRand, rhs.chainRand); 
  swap(proposals, rhs.proposals); 
}


void Chain::step()
{
#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*traln);
#endif

  currentGeneration++; 
  tree *tr = traln->getTr();   
  evaluator->imprint(*traln);

  // inform the rng that we produce random numbers for generation x  
  chainRand.rebase(currentGeneration);

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = getChainHeat();

  auto& pfun = drawProposalFunction();
 
  /* reset proposal ratio  */
  hastings = 0; 

#ifdef DEBUG_ACCEPTANCE
  double oldPrior = prior.getLnPrior();
#endif

  pfun->applyToState(*traln, prior, hastings, chainRand);
  pfun->evaluateProposal(evaluator, *traln, prior);
  
  double priorRatio = prior.getLnPriorRatio();
  double lnlRatio = tr->likelihood - prevLnl; 

  double testr = chainRand.drawRandDouble01();
  double acceptance = exp(( priorRatio + lnlRatio) * myHeat + hastings) ; 

#ifdef DEBUG_ACCEPTANCE
  tout  << endl << "(" << oldPrior<<  " - " << prior.getLnPriorRatio()  << ") + ( " << tr->likelihood << " - "<< prevLnl   << ") * " << myHeat << " + " << hastings << endl; 
#endif

  // hastings = 0; 

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

  evaluator->expensiveVerify(*traln); 
  
#ifdef DEBUG_TREE_LENGTH  
  assert( fabs (traln->getTreeLengthExpensive() - traln->getTreeLength())  < 1e-6); 
#endif

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*traln);
#endif

  // debug_checkTreeConsistency(this->getTraln());

  if(this->tuneFrequency <  pfun->getNumCallSinceTuning() )
    {
      pfun->autotune();
    }
}
