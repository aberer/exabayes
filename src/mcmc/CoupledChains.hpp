/**
 * @file CoupledChains.hpp
 *
 * represents a run (consisting of a number of coupled chains)
 *
 */


#ifndef _COUPLED_CHAINS_H
#define _COUPLED_CHAINS_H

#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "BlockRunParameters.hpp"
#include "Chain.hpp"
#include "SuccessCounter.hpp"
#include "TopologyFile.hpp"
#include "ParameterFile.hpp"
#include "SwapMatrix.hpp"
#include "SwapElem.hpp"

#include <queue>

class PendingSwap;
class ParallelSetup;

/**
 * @brief represents some coupled chains, one of them cold, many of
 * them hot
 */


///////////////////////////////////////////////////////////////////////////////
//                               COUPLED CHAINS                              //
///////////////////////////////////////////////////////////////////////////////
class CoupledChains : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    CoupledChains(
        Randomness   rand,
        int          runNum,
        string       workingdir,
        std::string  runname,
        int          numCoupled,
        vector<Chain>chains);
    // ________________________________________________________________________
    CoupledChains(
        CoupledChains&& rhs) = default;
    // ________________________________________________________________________
    CoupledChains(
        const CoupledChains& rhs) = default;
    // ________________________________________________________________________
    CoupledChains&                              operator=(
        const CoupledChains& rhs) = default;
    CoupledChains&                              operator=(
        CoupledChains&& rhs) = default;
    // ________________________________________________________________________
    void                                        run(
        uint64_t numGen);
    // ________________________________________________________________________
    void                                        executePart(
        uint64_t       startGen,
        uint64_t       numGen,
        ParallelSetup& pl);
    // ________________________________________________________________________
    void                                        doStep(
        nat           id,
        ParallelSetup&pl);
    // ________________________________________________________________________
    void                                        setSamplingFreq(
        nat i){_samplingFreq = i; }
    // ________________________________________________________________________
    void                                        setHeatIncrement(
        double temp){_heatIncrement = temp; }
    // ________________________________________________________________________
    void                                        setTemperature(
        double temp){_heatIncrement = temp;  }
    // ________________________________________________________________________
    std::vector<Chain>&                         getChains()
    {return _chains; }
    // ________________________________________________________________________
    nat                                         getRunid()  const
    {return _runid; }
    // ________________________________________________________________________
    const vector<Chain>&                        getChains() const
    {return _chains; }
    // ________________________________________________________________________
    size_t                                      getNumberOfChains()
    {return _chains.size(); }
    // ________________________________________________________________________
    void                                        setNumSwapsPerGen(
        double s){_numSwapsPerGen = s; }
    // ________________________________________________________________________
    void                                        setRunName(
        string a){_runname = a;  }
    // ________________________________________________________________________
    void                                        initializeOutputFiles(
        bool isDryRun);
    // ________________________________________________________________________
    SwapMatrix                                  getSwapInfo() const
    {return _swapInfo; }
    // ________________________________________________________________________
    void                                        setSwapInfo(
        SwapMatrix swap){_swapInfo = swap; }
    // ________________________________________________________________________
    void                                        addToSwapMatrix(
        const SwapMatrix&toAdd){_swapInfo = _swapInfo + toAdd;  }
    // ________________________________________________________________________
    const Randomness&                           getRandomness() const
    {return _rand; }
    // ________________________________________________________________________
    std::vector<std::string>                    getAllFileNames() const;
    // ________________________________________________________________________
    virtual void                                deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    void                                        regenerateOutputFiles(
        std::string _workdir,
        std::string prevId);
    // ________________________________________________________________________
    std::list<SwapElem>                         generateSwapsForBatch(
        uint64_t startGen,
        uint64_t numGen);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    bool                                        allMyChainsAreBlocked(
        const std::vector<bool>&isBlocked,
        const ParallelSetup&    pl) const;
    // ________________________________________________________________________
    bool                                        doLocalSwap(
        ParallelSetup& pl,
        const SwapElem&theSwap);
    // ________________________________________________________________________
    bool                                        doSwap(
        ParallelSetup&  pl,
        const SwapElem& elem);
    // ________________________________________________________________________
    PendingSwap                                 prepareSwap(
        ParallelSetup&  pl,
        const SwapElem& theSwap);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::vector<Chain>                    _chains;
    SwapMatrix                            _swapInfo;
    double                                _heatIncrement;
    Randomness                            _rand;
    int                                   _runid;
    int                                   _samplingFreq;
    string                                _runname;
    string                                _workdir;

    // order is coupled to the heat id
    std::unordered_map<nat, TopologyFile> _paramId2TopFile;
    std::vector<ParameterFile>            _pFile;

    double                                _numSwapsPerGen;
};


#endif

