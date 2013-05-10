#include <sstream>

#include "Chain.hpp"
#include "Topology.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "output.h"
#include "GlobalVariables.hpp"
#include "BipartitionHash.hpp"
#include "adapters.h"
#include "proposals.h"
#include "ExtendedTBR.hpp"
#include "ExtendedSPR.hpp"
#include "WrappedProposal.hpp"
#include "ParsimonySPR.hpp" 
#include "eval.h"
#include "TreeLengthMultiplier.hpp"
#include "StatNNI.hpp"
#include "PartitionProposal.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"


// meh =/ 
extern bool isNewProposal[NUM_PROPOSALS]; 

#include <sstream>

Chain::Chain(randKey_t seed, int id, int _runid, TreeAln* _traln, PriorBelief _prior, vector<Category> propCats) 
  : traln(_traln)
  , couplingId(id)
  , currentGeneration(0)
  , hastings(1)
  , proposalCategories(propCats)    
  , runid(_runid)
  , prior(_prior)
{
  chainRand = new Randomness(seed.v[0]);

  initParamDump(); 

  for(int j = 0; j < traln->getNumberOfPartitions(); ++j)
    traln->initRevMat(j);

  evaluateFullNoBackup(this);   
  
  tree *tr = traln->getTr();
  tout <<  "init lnl=" << traln->getTr()->likelihood << "\tTL=" << branchLengthToReal(tr, traln->getTreeLength()) << "\tseeds="  << *chainRand << endl; 
  saveTreeStateToChain(); 
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


/**
   @brief assigns the entire chain state from rhs to this chain

   This includes topology, parameters, all proposal parameters 
 */ 
Chain& Chain::operator=(Chain& rhs)
{
  TreeAln &thisTraln = *traln; 
  TreeAln &rhsTraln = *(rhs.traln); 

  thisTraln = rhsTraln; 

  // TODO should also copy proposals, but that is currently not used

  dump.topology->saveTopology(thisTraln);  
  return *this; 
}


void Chain::saveTreeStateToChain()
{
  dump.topology->saveTopology(*(traln));

  /* save model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = dump.infoPerPart + i; 
      pInfo *partition = traln->getPartition(i); 
      
      info->alpha =  partition->alpha; 
      
      memcpy(info->substRates, partition->substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, partition->frequencies, 4 * sizeof(double)); 
    }  
}




void Chain::applyChainStateToTree()
{
  /* TODO enable multi-branch    */
  assert(traln->getNumBranches() == 1); 

  dump.topology->restoreTopology(*(traln));

  /* restore model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = dump.infoPerPart + i ; 

      pInfo *partition = traln->getPartition( i);
      partition->alpha = info->alpha;

      memcpy(partition->substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(partition->frequencies, info->frequencies, 4 * sizeof(double));
      
      traln->initRevMat(i);
      traln->discretizeGamma(i);	 
    } 

  evaluateFullNoBackup(this); 
}




void Chain::debug_printAccRejc(AbstractProposal *prob, bool accepted, double lnl) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  if(isOutputProcess())
    {
      if(accepted)
	printInfo( "ACC\t");   
      else 
	printInfo("rej\t");   	  
      printf("%s %g\n" ,prob->getName().c_str() , lnl); 
    }
#endif
}




AbstractProposal* Chain::drawProposalFunction()
{    
  double r = chainRand->drawRandDouble01();  
  for(auto c : proposalCategories)
    {
      double ref = c.getCatFreq(); 
      if(r <= ref )
	return c.drawProposal(*chainRand);
      else 
	r -= ref; 
    }
  
  assert(0); 
  return NULL; 
}



void Chain::step()
{
  this->currentGeneration++; 

  this->restorer->resetRestorer(*traln);

  // inform the rng that we produce random numbers for generation x  
  chainRand->rebase(currentGeneration);

  tree *tr = this->traln->getTr();   

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = getChainHeat();

  AbstractProposal *pfun = drawProposalFunction();
 
  /* reset proposal ratio  */
  this->hastings = 1; 
  
  double oldPrior = prior.getLogProb();
  pfun->applyToState(*traln, prior, hastings, *chainRand);
  pfun->evaluateProposal(*traln, prior);
  
  double priorRatio = prior.getLogProb() - oldPrior;
  double lnlRatio = tr->likelihood - prevLnl; 

  double testr = chainRand->drawRandDouble01();
  double acceptance = 
    exp(( priorRatio   + lnlRatio) * myHeat) 
    * hastings;
  
  bool wasAccepted  = testr < acceptance; 

  // printInfo("prior=%f\tlnl=%f\tacc=%f\ttestr=%f\n", priorRatio, lnlRatio, acceptance, testr); 
  
  debug_printAccRejc( pfun, wasAccepted, tr->likelihood); 


#ifdef VERIFY_LNL_SUPER_EXPENSIVE  
  // TEST
  double lnlBefore = tr->likelihood; 
  evaluateFullNoBackup(this); 
  assert( ( this->traln->getTr()->likelihood  - lnlBefore ) < ACCEPTED_LIKELIHOOD_EPS); 
  // END
#endif


  if(wasAccepted)
    {
      pfun->accept();      
      expensiveVerify(this);
    }
  else
    {
      pfun->resetState(*traln, prior);
      pfun->reject();

      // TODO maybe not here =/ 
      this->restorer->restoreArrays(*traln); // restores the previous tree state 
    }

#ifdef VERIFY_LNL_SUPER_EXPENSIVE  
  // {
  //   double lnlBefore = tr->likelihood;  
  //   evaluateFullNoBackup(this);
  //   assert( ( this->traln->getTr()->likelihood  - lnlBefore ) < ACCEPTED_LIKELIHOOD_EPS); 
  // }
#endif

  debug_checkTreeConsistency(this->traln->getTr());
 
  if(this->tuneFrequency <  pfun->getNumCallSinceTuning() )
    pfun->autotune();


#if 0 
  /* autotuning for proposal parameters. With increased parallelism
     this will become more complicated.  */
  if( globals.tuneFreq > 0 && this->currentGeneration % globals.tuneFreq == globals.tuneFreq - 1 )
    {
      for(int i = 0; i < this->numProposals; ++i)
	{
	  proposalFunction *pf = this->proposals[i]; 
	  
	  if(pf->autotune)	/* only, if we set this   */
	    {
#ifdef TUNE_ONLY_IF_ENOUGH
	      if(pf->sCtr.lAcc + pf->sCtr.lRej < globals.tuneFreq )
		continue; 
#endif
	      pf->autotune(this, pf);
	    }
	}
    }
#endif


}





void Chain::initParamDump()
{  
  tree *tr = traln->getTr();
  
  dump.topology = new Topology(tr->mxtips); 
  dump.infoPerPart = (perPartitionInfo*)exa_calloc(traln->getNumberOfPartitions(), sizeof(perPartitionInfo));
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *p =  dump.infoPerPart + i ; 
      p->alpha = 0.5; 
      for(int j = 0; j < 6; ++j)
	p->substRates[j] = 0.5; 
      p->substRates[5] = 1; 
      for(int j = 0; j < 4; ++j)
	p->frequencies[j] = 0.25;       
    }
}


/**
   @brief prints a message with associated chain/run/heat information. 

   Please always use this, when possible, it also assures that the
   message is only printed once.
 */ 
// void Chain::printInfo(const char *format, ...)
// {  
//   if( isOutputProcess())
//     {
//       printf("[run %d / heat %d / gen %d] ", this->runid, this->couplingId, this->currentGeneration); 
//       va_list args;
//       va_start(args, format);     
//       vprintf(format, args );
//       va_end(args);
//     }
// }




void Chain::clarifyOwnership()
{
  for(auto c : proposalCategories)
    for(auto p : c.getProposals())
      p->setOwningChain(this); 
}


void Chain::printParams(FILE *fh)
{
  tree *tr = traln->getTr(); 


  double treeLength = branchLengthToReal(tr, traln->getTreeLength()); 
  assert(treeLength != 0.); 
  fprintf(fh, "%d\t%f\t%.3f", currentGeneration,
	  tr->likelihood,  
	  treeLength); 

  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      pInfo *partition = traln->getPartition(i); 
      fprintf(fh, "\t%f\t%f", traln->accessPartitionLH( i),traln->getAlpha(i)) ; 
      for(int j = 0; j < 6 ; ++j) /* TODO */
	fprintf(fh, "\t%.2f", partition->substRates[j]); 
      for(int j = 0; j < 4 ; ++j) /* TODO */
	fprintf(fh, "\t%.2f", traln->getFrequency(i,j));
    }

  fprintf(fh, "\n"); 
  fflush(fh); 
}




void Chain::printParamFileStart(FILE *fh)
{
  char *tmp = "TODO"; 
  fprintf(fh, "[ID: %s]\n", tmp);
  fprintf(fh, "Gen\tLnL\tTL"); 

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
  tree *tr = traln->getTr();
  memset(traln->getTr()->tree_string, 0, traln->getTr()->treeStringLength * sizeof(char) ); 
  
  Tree2stringNexus(traln->getTr()->tree_string, tr,  traln->getTr()->start->back, 0 ); 
  fprintf(fh,"\ttree gen.%d = [&U] %s\n", currentGeneration, traln->getTr()->tree_string);
  fflush(fh);
}


void Chain::printNexusTreeFileStart( FILE *fh  )
{
  fprintf(fh, "#NEXUS\n[ID: %s]\n[Param: tree]\nbegin trees;\n\ttranslate\n", "TODO" ); 

  for(int i = 0; i < traln->getTr()->mxtips-1; ++i)
    fprintf(fh, "\t\t%d %s,\n", i+1, traln->getTr()->nameList[i+1]); 
  fprintf(fh, "\t\t%d %s;\n", traln->getTr()->mxtips, traln->getTr()->nameList[traln->getTr()->mxtips]);
}
