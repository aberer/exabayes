#include "ProposalFactory.hpp"

ProposalFactory::ProposalFactory()
  : branchLengthMultiplier ( 1.386294)
  , rateSlidingWindow ( 0.15)
  , frequencySlidingWindow ( 0.2 )
  , gammaSlidingWindow ( 0.75)
  , secondaryBranchLengthMultiplier ( 0.098)
  , treeLengthMultiplier ( 1.386294)
  , dirichletAlpha ( 100 )
  , gammaMultiplier ( 0.811 )
  , nodeSliderMultiplier ( 0.191 )
{
}



vector<ProposalPtr> ProposalFactory::getProposals(const vector<RandomVariablePtr> &variables) const
{
  assert(0); 
  vector<ProposalPtr> result; 
  return result; 
}
