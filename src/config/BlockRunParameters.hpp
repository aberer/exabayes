#ifndef _BLOCK_RUNPARAMETERS_H
#define _BLOCK_RUNPARAMETERS_H

#include "ExaBlock.hpp"
#include "common.h"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////
//                            BLOCK RUN PARAMETERS                           //
///////////////////////////////////////////////////////////////////////////////
class BlockRunParameters : public ExaBlock
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BlockRunParameters();
    // ________________________________________________________________________
    BlockRunParameters(
        const BlockRunParameters& rhs) = default;
    // ________________________________________________________________________
    BlockRunParameters(
        BlockRunParameters&& rhs) = default;
    // ________________________________________________________________________
    BlockRunParameters&                   operator=(
        const BlockRunParameters&rhs) = default;
    // ________________________________________________________________________
    BlockRunParameters&                   operator=(
        BlockRunParameters&&rhs) = default;
    // ________________________________________________________________________
    virtual void                          Read(
        NxsToken&token);
    // ________________________________________________________________________
    nat                                   getTuneFreq() const
    {return tuneFreq; }
    // ________________________________________________________________________
    bool                                  getTuneHeat() const
    {return tuneHeat; }
    // ________________________________________________________________________
    double                                getNumSwapsPerGen() const
    {return numSwapsPerGen; }
    // ________________________________________________________________________
    double                                getHeatFactor() const
    {return heatFactor; }
    // ________________________________________________________________________
    nat                                   getPrintFreq() const
    {return printFreq; }
    // ________________________________________________________________________
    nat                                   getNumCoupledChains() const
    {return numCoupledChains; }
    // ________________________________________________________________________
    bool                                  isUseAsdsfMax() const
    {return useAsdsfMax; }
    // ________________________________________________________________________
    uint64_t                              getNumGen() const
    {return numGen; }
    // ________________________________________________________________________
    nat                                   getNumRunConv() const
    {return numRunConv; }
    // ________________________________________________________________________
    nat                                   getSamplingFreq() const
    {return samplingFreq; }
    // ________________________________________________________________________
    double                                getBurninProportion() const
    {return burninProportion; }
    // ________________________________________________________________________
    nat                                   getBurninGen() const
    {return burninGen; }
    // ________________________________________________________________________
    double                                getAsdsfIgnoreFreq() const
    {return asdsfIgnoreFreq; }
    // ________________________________________________________________________
    nat                                   getDiagFreq() const
    {return diagFreq; }
    // ________________________________________________________________________
    double                                getAsdsfConvergence() const
    {return asdsfConvergence; }
    // ________________________________________________________________________
    bool                                  isUseParsimonyStarting() const
    {return useParsimonyStarting; }
    // ________________________________________________________________________
    bool                                  isHeatedChainsUseSame() const
    {return heatedChainsUseSame; }
    // ________________________________________________________________________
    nat                                   getChkpntFreq() const
    {return chkpntFreq; }
    // ________________________________________________________________________
    bool                                  isComponentWiseMH() const
    {return componentWiseMH; }
    // ________________________________________________________________________
    bool                                  isUseStopCriterion() const
    {return useStopCriterion; }
    // ________________________________________________________________________
    void                                  verify() const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    nat      diagFreq;
    double   asdsfIgnoreFreq;
    double   asdsfConvergence;
    bool     useStopCriterion;
    nat      burninGen;
    double   burninProportion;
    int      samplingFreq;
    int      numRunConv;
    uint64_t numGen;
    int      numCoupledChains;
    int      printFreq;
    double   heatFactor;
    bool     tuneHeat;
    nat      tuneFreq;
    bool     useParsimonyStarting;
    bool     heatedChainsUseSame;
    nat      chkpntFreq;
    bool     componentWiseMH;
    bool     useAsdsfMax;
    double   numSwapsPerGen;
};


#endif
