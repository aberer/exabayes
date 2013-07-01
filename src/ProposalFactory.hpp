#ifndef _PROPOSAL_FACTORY_H
#define _PROPOSAL_FACTORY_H

#include "axml.h"
#include "AbstractProposal.hpp"

class ProposalFactory
{
public: 
  ProposalFactory(); 
  vector<ProposalPtr> getProposals(const vector<RandomVariablePtr> &variables) const; 
  
private: 			// ATTRIBUTES
  double frequencySlidingWindow ; 
  double branchLengthMultiplier; 
  double rateSlidingWindow; 
  double gammaSlidingWindow; 
  double secondaryBranchLengthMultiplier; 
  double dirichletAlpha; 
  double treeLengthMultiplier ; 
  double gammaMultiplier; 
  double nodeSliderMultiplier; 
  

}; 


#endif
