#include <sstream>


#include "Chain.hpp"
#include "Topology.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "output.h"
#include "GlobalVariables.hpp"
#include "BipartitionHash.hpp"
#include "tune.h"
#include "ExtendedTBR.hpp"
#include "ExtendedSPR.hpp"
#include "ParsimonySPR.hpp" 
#include "eval.h"
#include "TreeLengthMultiplier.hpp"
#include "StatNNI.hpp"
#include "PartitionProposal.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"


// #define DEBUG_ACCEPTANCE

Chain::Chain(randKey_t seed, int id, int _runid, shared_ptr<TreeAln> _traln, const vector< unique_ptr<AbstractProposal> > &_proposals, int _tuneFreq ,const vector<RandomVariable> &variables) 
  : traln(_traln)
  , runid(_runid)
  , tuneFrequency(_tuneFreq)
  , hastings(0)
  , currentGeneration(0)
  , couplingId(id)
  , state(*traln)
  , chainRand(seed.v[0])
  , relWeightSum(0)
  , prior(*_traln, variables)
{
  for(auto &p : _proposals)
    proposals.push_back(unique_ptr<AbstractProposal>(p->clone())); 

  for(int j = 0; j < traln->getNumberOfPartitions(); ++j)
    traln->initRevMat(j);

  // evaluateFullNoBackup(*traln);   

  // addChainInfo(tout)  << " lnPr="  << prior.getLnPrior() << " lnLH=" << traln->getTr()->likelihood << "\tTL=" << branchLengthToReal(traln->getTr(), traln->getTreeLengthExpensive()) << "\tseeds=>"  << chainRand << endl; 

  saveTreeStateToChain(); 

  for(auto& elem : proposals)
    relWeightSum +=  elem->getRelativeWeight();
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

  evaluateFullNoBackup(*traln); 
  prior.reinitPrior(*traln);
}




void Chain::debug_printAccRejc(unique_ptr<AbstractProposal> &prob, bool accepted, double lnl, double lnPr ) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  if(isOutputProcess())
    tout << "[run=" << runid << ",heat="  << couplingId << ",gen="  << currentGeneration << "]\t" << (accepted ? "ACC" : "rej" )  << "\t"<< prob->getName() << "\t" << setprecision(2) << fixed << lnl  <<  "\t"<< lnPr << endl; 
#endif
}




unique_ptr<AbstractProposal>& Chain::drawProposalFunction()
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

  double treeLength =  branchLengthToReal(tr, traln->getTreeLengthExpensive()); 
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
  swap(chainRand, rhs.chainRand); 
  swap(proposals, rhs.proposals); 
}


void Chain::step()
{
  // cout << "current prior is " << prior.getLnPrior() << " and ratio is " << prior.getLnPriorRatio() << endl; 
  prior.verifyPrior(*traln);

  currentGeneration++; 
  tree *tr = traln->getTr();   
  traln->getRestorer()->resetRestorer(*traln);

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
  pfun->evaluateProposal(*traln, prior);
  
  double priorRatio = prior.getLnPriorRatio();
  double lnlRatio = tr->likelihood - prevLnl; 

  double testr = chainRand.drawRandDouble01();
  double acceptance = exp(( priorRatio   + lnlRatio) * myHeat + hastings) ; 

#ifdef DEBUG_ACCEPTANCE
  tout << "(" << prior.getLogProb() << " - " << oldPrior << ") + ( " << tr->likelihood << " - "<< prevLnl   << ") * " << myHeat << " + " << hastings << endl; 
#endif

  bool wasAccepted  = testr < acceptance; 

  debug_printAccRejc( pfun, wasAccepted, tr->likelihood, prior.getLnPrior()); 

  if(wasAccepted)
    {
      pfun->accept();      
      prior.accept();
    }
  else
    {
      pfun->resetState(*traln, prior);
      pfun->reject();
      prior.reject();

      traln->getRestorer()->restoreArrays(*traln);
    }

  expensiveVerify(*traln); 

// proposals/Chai
  // TreePrinter tp(false, true, false); 
  // cout << "after proposal: " << tp.printTree(*traln)<< endl; 

  

#ifdef DEBUG_TREE_LENGTH  
  assert( fabs (traln->getTreeLengthExpensive() - traln->getTreeLength())  < 1e-6); 
#endif

#ifdef DEBUG_VERIFY_LNPR
  prior.verifyPrior(*traln);
#endif

  debug_checkTreeConsistency(this->getTraln());
 
  if(this->tuneFrequency <  pfun->getNumCallSinceTuning() )
    pfun->autotune();

}
