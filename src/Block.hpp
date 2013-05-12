#ifndef _BLOCK_H 
#define _BLOCK_H 

#include <map>

#include <ncl/ncl.h>

#include "SampleMaster.hpp"

using namespace std; 




class ExabayesBlock : public NxsBlock
{
public: 
  ExabayesBlock(SampleMaster *_sm ) ; 

  virtual void Read(NxsToken &token); 
  void setupMap(); 
  void fillProposalWeights(vector<double> &weights){ weights = proposalWeights;  } 
  PriorBelief getPrior() {return prior; }

private: 
  SampleMaster *sm; 
  vector<double> proposalWeights ; 
  PriorBelief prior;   
  map<string, proposal_type> name2proposal; 

  shared_ptr<AbstractPrior> parsePrior(NxsToken &token, NxsString &value);
}; 


#endif
