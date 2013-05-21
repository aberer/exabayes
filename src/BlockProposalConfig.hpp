#ifndef _BLOCK_PROPOSALCONFIG_H
#define _BLOCK_PROPOSALCONFIG_H

#include <ncl/ncl.h>
#include <map>

#include "GlobalVariables.hpp"

class BlockProposalConfig : public NxsBlock
{
public: 
  BlockProposalConfig();
  virtual void Read(NxsToken &token); 
  void setupMap(); 
  void fillProposalWeights(vector<double> &weights){ weights = proposalWeights;  } 

  double getEsprStopProp() const {return esprStopProp; } 
  double getParsimonyWarp() const {return parsimonyWarp; }
  int getGuidedRadius() const {return guidedRadius; }

private: 
  vector<double> proposalWeights ; 
  map<string, proposal_type> name2proposal; 

  double esprStopProp; 
  double parsimonyWarp;   
  int guidedRadius; 



}; 


#endif
