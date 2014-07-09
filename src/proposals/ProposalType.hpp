#ifndef _PROPOSAL_TYPE 
#define _PROPOSAL_TYPE 

#include <unordered_map>
#include <vector>

#include "model/Category.hpp" 


// add your new proposal here and update accordingly in everything in
// ProposalType.cpp .  if you missed something, then exabayes usually
// complains with an > assert(0);

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
    BRANCH_DIST_GAMMA = 17,
    DIRICH_REVMAT_ALL = 18, 
    LIKE_SPR = 19 ,
    DIRICH_REVMAT_PER_RATE = 20,
    SLIDING_REVMAT_PER_RATE = 21, 
    BL_DIST_WEIBULL = 22, 
    BL_SLID_GAMMA = 23 , 
    DIV_TIME_DIRICH = 24 
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
  /** 
      @brief gets the verbose name of the proposal 
   */ 
  std::string getLongName(ProposalType type); 
  /** 
      @brief gets the proposal type from the name in the config file   
   */ 
  ProposalType getTypeFromConfigString(std::string s); 
  /** 
      @brief gets the name of the proposal by type 
   */ 
  std::string getConfigStringFromType(ProposalType p ); 
  /** 
      @brief gets all proposals for a category that integrate only over one parameter (primarily) 
   */ 
  std::vector<ProposalType> getSingleParameterProposalsForCategory(Category c) ; 
  /** 
      @brief IMPORTANT get all relevant proposals for a category  
   */ 
  // std::vector<ProposalType> getProposalsForCategory(Category c) ; 
  /** 
      @brief gets all relevant proposals 
   */ 
  std::vector<ProposalType> getAllProposals();   
  /**
     @brief indicates, if the string specifies a valid proposal name 
   */ 
  bool isValidName(std::string name); 
  /** 
      @brief indicates, if the proposal is ready for productive use
      (if this is not the case, then the proposal can only be
      activated by explicitly specifying it in the config file)
   */ 
  bool isReadyForProductiveUse(ProposalType p); 

}

#endif
