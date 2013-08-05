#ifndef _BLOCK_RUNPARAMETERS_H
#define _BLOCK_RUNPARAMETERS_H


#include <cassert>
#include <ncl/ncl.h>

#include "common.h"


class BlockRunParameters : public NxsBlock
{
public: 
  BlockRunParameters(); 

  virtual void Read(NxsToken &token); 

  // getters 
  nat getTuneFreq() const { return tuneFreq;  }
  bool getTuneHeat() const { return tuneHeat; }
  nat getSwapInterval() const { return swapInterval; }
  double getHeatFactor() const { return heatFactor ; }
  nat getPrintFreq() const { return printFreq; }
  nat getNumCoupledChains() const { return numCoupledChains; }
  // string getRunId() const { return runId; }
  nat getNumGen() const { return numGen; }
  nat getNumRunConv() const { return numRunConv; }
  nat getSamplingFreq() const { return samplingFreq; }
  double getBurninProportion() const { return burninProportion; }
  nat getBurninGen() const { return burninGen; }
  double getAsdsfIgnoreFreq() const { return asdsfIgnoreFreq; 	}
  nat getDiagFreq() const { return diagFreq ; }
  double getAsdsfConvergence() const {return asdsfConvergence; }
  bool isUseParsimonyStarting() const {return useParsimonyStarting; } 
  bool isHeatedChainsUseSame() const {return heatedChainsUseSame;}
  nat getChkpntFreq() const {return chkpntFreq; }
  bool isComponentWiseMH() const {return componentWiseMH; }

private: 
  nat diagFreq ; 
  double asdsfIgnoreFreq; 	
  double asdsfConvergence; 
  nat burninGen; 
  double burninProportion; 
  nat samplingFreq; 
  nat numRunConv; 
  nat numGen; 
  nat numCoupledChains; 
  nat printFreq; 
  double heatFactor ; 
  nat swapInterval; 
  bool tuneHeat; 
  nat tuneFreq;  
  bool useParsimonyStarting; 
  bool heatedChainsUseSame; 
  nat chkpntFreq; 
  bool componentWiseMH; 
}; 


#endif
