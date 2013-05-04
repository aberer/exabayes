#ifndef _PROPOSAL_TYPE_H 
#define  _PROPOSAL_TYPE_H 


/* okay, so defining enums this way is rather save  */
#define NUM_PROPOSALS (23) //PROPOSALADD NUM_PROPOSALS NOTE Do not remove/modify  this line except for numerical value. The script addProposal.pl needs it as an identifier.
typedef enum
  {
    TL_MULT = 0, 
    UPDATE_MODEL = 1 ,
    UPDATE_GAMMA = 2 ,
    UPDATE_GAMMA_EXP=3,
    UPDATE_SINGLE_BL = 4,
    UPDATE_SINGLE_BL_EXP = 5 ,
    UPDATE_SINGLE_BL_BIUNIF = 6,
    UPDATE_MODEL_BIUNIF = 7,
    UPDATE_MODEL_SINGLE_BIUNIF = 8,
    UPDATE_MODEL_ALL_BIUNIF = 9,
    UPDATE_MODEL_PERM_BIUNIF = 10,
    UPDATE_FREQUENCIES_BIUNIF = 11,
    BRANCH_LENGTHS_MULTIPLIER =12,
    FREQUENCY_SLIDER = 13,
    GUIDED_SPR = 14,
    ST_NNI= 15, 
    GAMMA_MULTI = 16,
    NODE_SLIDER = 17,
    UPDATE_FREQUENCIES_DIRICHLET = 18,
    UPDATE_MODEL_DIRICHLET = 19,
    UPDATE_SINGLE_BL_GUIDED = 20,
    E_SPR = 21,
    E_TBR = 22
    //PROPOSALADD proposal_type NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.


  } proposal_type;


#endif
