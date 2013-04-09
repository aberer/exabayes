#ifndef _PROPOSAL_TYPE_H 
#define  _PROPOSAL_TYPE_H 


/* okay, so defining enums this way is rather save  */
#define NUM_PROPOSALS (17) //PROPOSALADD NUM_PROPOSALS NOTE Do not remove/modify  this line except for numerical value. The script addProposal.pl needs it as an identifier.
typedef enum
{
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
E_SPR_MAPPED = 12,
BRANCH_LENGTHS_MULTIPLIER =13,
FREQUENCY_SLIDER = 14,
GUIDED_SPR = 15,
ST_NNI= 16
//PROPOSALADD proposal_type NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.


} proposal_type;

#endif
