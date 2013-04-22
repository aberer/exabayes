#ifndef _PRSF_H
#define _PRSF_H


void printPRSF(char *runId ); 
int guessNumGen( int numChain, int numParam, char *runid); 
int guessNumChain(char *runid);

int guessNumParam(char *runid); 
void parseFiles(int nGens, int mChains, int numParam, char *runId, double ****resultPtr, char ***paramNamesPtr); 
double getPrsfForParameter(int paramNum, int numChain, int numGen, double ***matrix); 
void printMatrix(int numParam, int nGens, int mChains, char **paramNames, double ***result ); 
void freeMatrix(int mChains, int numParam, double ***result); 
#endif
