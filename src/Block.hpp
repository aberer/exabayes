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

private: 
  SampleMaster *sm; 
  vector<double> proposalWeights ; 
  map<string, proposal_type> name2proposal; 
}; 


#endif
