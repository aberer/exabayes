#pragma once 

#include "axml.h"
#include "bayes.h"

typedef node* nodep;  
typedef struct _state state; 



class LnlRestorer
{
public: 
  LnlRestorer(state *chain);
  void resetRestorer(); 
  void restore();  
  void traverseAndSwitchIfNecessary(nodep virtualRoot, int model, bool fullTraversal); 
  ~LnlRestorer();

private:   
  void storeOrientation(); 
  void loadOrientation();
  void swapArray(int number, int model); 

  int modelEvaluated; 
  state *chain; 
  double ***reserveArrays; 
  int* orientation; 
  bool *wasSwitched; 
  nat **partitionScaler;   
  double prevLnl; 		// for DEBUG
}; 



