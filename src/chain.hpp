#pragma once 

class TreeAln; 
class LnlRestorer; 
class Topology; 


typedef struct
{
  int modelNum; 
  double alpha ;
 
  /* TODO only works with DNA */
  int numRates; 
  double substRates[6]; 

  int numFreqs; 
  double frequencies[4];

} perPartitionInfo; 		/* relevant info from pInfo  */


typedef struct 
{
  /* topology */
  // topol *topo; 
  Topology *topology; 
  
  /* branch lengths  */
  double *branchLengths; 

  /* substitution parameter */
  perPartitionInfo *infoPerPart; 

} paramDump; 


typedef struct _pfun proposalFunction; 


typedef struct _state
{  
  TreeAln *traln; 
  /* tree* tr; */

  boolean wasAccepted; 	/// for debug only  

  int id;   
  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  int currentGeneration;   
  
  double priorProb; /// the prior probability of the current state   => store in log form  

  LnlRestorer *restorer; 

  /* lnlContainer lnl;  /// the current likelihood of the state */

  double hastings;/// the proposal ratio 

  double penaltyFactor; //change to the probability of picking a proposal  

  proposalFunction **proposals; 
  int numProposals; 
  double *categoryWeights; 

  /* new stuff that we need when having multiple chains  */
  FILE *topologyFile; 
  FILE *outputParamFile; 

  /* RNG */
  randKey_t rKey;
  randCtr_t rCtr;

  /* saves the entire space in the parameter space  */
  paramDump dump;
} state;

