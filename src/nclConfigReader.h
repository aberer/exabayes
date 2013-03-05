
#ifndef _NCL_CONFIG_READER
#define _NCL_CONFIG_READER


typedef struct  
{
  double initSPRWeight;
  double initGammaWeight; 
  double initGammaExpWeight;
  double initModelWeight;
  double initSingleBranchWeight; 
  double initSingleBranchExpWeight;
double initSingleBranchBiunifWeight;
double initModelBiunifWeight;
double initModelSingleBiunifWeight;
double initModelAllBiunifWeight;
double initModelPermBiunifWeight;
double initFrequenciesWeight;
double initEsprMappedWeight;
  //PROPOSALADD initParamStruct NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
  double initPenaltyFactor;   
  double eSprStopProb; 
  int numGen; 
  int samplingFrequency; 
  int numIndiChains; 
  int diagFreq; 

} initParamStruct ; 


#ifdef __cplusplus
extern "C"
{
#endif
  void parseConfigWithNcl(char *configFileName, initParamStruct **initParam); 
#ifdef __cplusplus
}
#endif

#endif
