#ifndef _PROPOSAL_TYPE 
#define _PROPOSAL_TYPE 

#include <unordered_map>
#include <vector>

#include "Category.hpp" 

enum class ProposalType
{
    ST_NNI= 0, 
    E_SPR = 1,
    E_TBR = 2,
    PARSIMONY_SPR= 3 , 
    GUIDED_SPR = 4,
    BRANCH_SLIDER = 5 ,
    TL_MULT = 6, 
    BRANCH_COLLAPSER = 7,
    NODE_SLIDER = 8,
    BRANCH_LENGTHS_MULTIPLIER = 9 , 
    REVMAT_SLIDER = 10 ,
    REVMAT_DIRICHLET = 11, 
    RATE_HET_SLIDER = 12, 
    RATE_HET_MULTI = 13,     
    FREQUENCY_SLIDER = 14, 
    FREQUENCY_DIRICHLET = 15,     
    AMINO_MODEL_JUMP = 16,
    BRANCH_GIBBS = 17  
}; 


class ProposalTypeHash
{
public: 
  size_t operator()(const ProposalType& t) const 
  {
    std::hash<int> hasher; 
    return hasher(int(t)) ; 
  }
}; 


namespace ProposalTypeFunc
{
  std::string getLongName(ProposalType type); 
  ProposalType getTypeFromConfigString(std::string s); 
  std::string getConfigStringFromType(ProposalType p ); 

  std::vector<ProposalType> getProposalsForCategory(Category c) ; 
  std::vector<ProposalType> getAllProposals(); 

  bool isValidName(std::string name); 
  bool isReadyForProductiveUse(ProposalType p); 

}; 

#endif
