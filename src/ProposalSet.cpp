#include "ProposalSet.hpp"


ProposalSet::ProposalSet(double relWeight, std::vector<std::unique_ptr<AbstractProposal> > _proposals)
  : relativeWeight(relWeight)
{
  for(auto &p : _proposals)
    proposals.emplace_back(p->clone()); 
  
#ifdef UNSURE
  assert(0);
#endif
  // TODO should not be necessary  
  for(auto &p : proposals)
    p->setInSetExecution(true);
} 


ProposalSet::ProposalSet(const ProposalSet &rhs)
  : relativeWeight(rhs.relativeWeight)
{
  for(auto &p : rhs.proposals)
    proposals.emplace_back(p->clone());

#ifdef UNSURE
  assert(0);
#endif
  // TODO should not be necessary  
  for(auto &p : rhs.proposals)
    p->setInSetExecution(true); 
}


ProposalSet& ProposalSet::operator=(ProposalSet rhs)
{
  std::swap(rhs, *this); 
  return *this; 
}


std::vector<AbstractProposal*> ProposalSet::getProposalView() const
{
  std::vector<AbstractProposal*> result; 
  for(auto &p : proposals)
    result.push_back(p.get());
  return result; 
}  


void ProposalSet::printVerboseAbbreviated(std::ostream &out, double sum)
{
  out << relativeWeight / sum * 100 <<   "%\tSET(totalNumber=" << proposals.size() << "):" << std::endl;
  for(auto &p : proposals)
    {
      // out << p.get() << std::endl << std::endl;   
      out << "\t"; 
      p->printShort(out); 
      out << std::endl; 
    }
}



void ProposalSet::readFromCheckpoint( std::istream &in )
{
  for(auto &p : proposals)
    p->readFromCheckpoint(in);    
}


void ProposalSet::writeToCheckpoint( std::ostream &out) const
{
  for(auto &p : proposals)
    p->writeToCheckpoint(out); 
}   


std::ostream& operator<<(std::ostream& out, const ProposalSet &rhs)
{
  out << "SET("; 
  
  std::unordered_set<std::string> pNames; 
  for(auto &p : rhs.proposals)
    pNames.insert(p->getName()); 

  out << pNames << ")"; 

  return out; 
}


bool ProposalSet::needsFullTraversal()
{
  bool result = true; 
  for(auto &p : proposals)
    result &= p->isNeedsFullTraversal();
  return result; 
}
