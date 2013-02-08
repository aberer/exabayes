
#ifndef _NCL_CONFIG_READER
#define _NCL_CONFIG_READER


typedef struct  
{
  double initSPRWeight;
  double initGammaWeight; 
  double initModelWeight;
  double initSingleBranchWeight; 
  int numGen; 
  double initPenaltyFactor; 
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
