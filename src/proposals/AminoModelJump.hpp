#ifndef _AAMODEL_JUMP_H  	// use whatever you want as an include guard (must NOT occur twice! )
#define _AAMODEL_JUMP_H

#include <vector>
using namespace std; 

#include "AbstractProposal.hpp"
#include "axml.h"

enum aaMatrix_t
  { 
    DAYHOFF_T =     DAYHOFF,
    DCMUT_T =      DCMUT,
    JTT_T =        JTT,
    MTREV_T =      MTREV,
    WAG_T =        WAG,
    RTREV_T =      RTREV,
    CPREV_T =      CPREV,
    VT_T =         VT,
    BLOSUM62_T =   BLOSUM62,
    MTMAM_T =      MTMAM,
    LG_T =         LG,
    MTART_T =      MTART,
    MTZOA_T =      MTZOA,
    PMB_T =        PMB,
    HIVB_T =       HIVB,
    HIVW_T =       HIVW,
    JTTDCMUT_T =   JTTDCMUT,
    FLU_T =        FLU
  }; 

class AminoModelJump : public AbstractProposal
{
public: 
  AminoModelJump(vector<aaMatrix_t> matrices); 

  virtual void readFromCheckpointCore(std::ifstream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ofstream &out)const { } //disabled

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior); 
  virtual void autotune()  ;
  virtual AbstractProposal* clone() const ;  

private: 
  vector<aaMatrix_t> matrices; 
}; 

#endif
