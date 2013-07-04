#ifndef _BLOCK_RUNPARAMETERS_H
#define _BLOCK_RUNPARAMETERS_H


#include <cassert>
#include <ncl/ncl.h>


class BlockRunParameters : public NxsBlock
{
public: 
  BlockRunParameters(); 

  virtual void Read(NxsToken &token); 

  // getters 
  int getTuneFreq() const { return tuneFreq;  }
  bool getTuneHeat() const { return tuneHeat; }
  int getSwapInterval() const { return swapInterval; }
  double getHeatFactor() const { return heatFactor ; }
  int getPrintFreq() const { return printFreq; }
  int getNumCoupledChains() const { return numCoupledChains; }
  string getRunId() const { return runId; }
  int getNumGen() const { return numGen; }
  int getNumRunConv() const { return numRunConv; }
  int getSamplingFreq() const { return samplingFreq; }
  double getBurninProportion() const { return burninProportion; }
  int getBurninGen() const { return burninGen; }
  double getAsdsfIgnoreFreq() const { return asdsfIgnoreFreq; 	}
  int getDiagFreq() const { return diagFreq ; }
  double getAsdsfConvergence() const {return asdsfConvergence; }
  bool isUseParsimonyStarting() const {return useParsimonyStarting; } 
  bool isHeatedChainsUseSame() const {return heatedChainsUseSame;}

private: 
  int diagFreq ; 
  double asdsfIgnoreFreq; 	
  double asdsfConvergence; 
  int burninGen; 
  double burninProportion; 
  int samplingFreq; 
  int numRunConv; 
  int numGen; 
  string runId; 
  int numCoupledChains; 
  int printFreq; 
  double heatFactor ; 
  int swapInterval; 
  bool tuneHeat; 
  int tuneFreq;  
  bool useParsimonyStarting; 
  bool heatedChainsUseSame; 
}; 


#endif
