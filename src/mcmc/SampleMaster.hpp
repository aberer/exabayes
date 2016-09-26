/**
 *  @file SampleMaster.hpp
 *
 *  @brief represents various chains sampling the posterior probability space
 *
 *  Despite of its modest name, this is in fact the master class.
 */

#ifndef _SAMPLING_H
#define _SAMPLING_H

#include "ProposalSet.hpp"
#include "BlockRunParameters.hpp"
#include "BlockProposalConfig.hpp"
#include "CommandLine.hpp"
#include "CoupledChains.hpp"
#include "ConfigReader.hpp"
#include "ParallelSetup.hpp"
#include "TimeTracker.hpp"
#include "Serializable.hpp"
#include "DiagnosticsFile.hpp"

#include <vector>

///////////////////////////////////////////////////////////////////////////////
//                               SAMPLE MASTER                               //
///////////////////////////////////////////////////////////////////////////////
class SampleMaster : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                          PUBLIC INTERFACE                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    SampleMaster(
        size_t numCpu);
    // ________________________________________________________________________
    SampleMaster(
        const SampleMaster& rhs) = default;
    // ________________________________________________________________________
    SampleMaster&
                        operator=(
        const SampleMaster&rhs)  = default;
    // ________________________________________________________________________
    /**
     *  @brief initializes the runs
     *  @notice this is the top-level function
     */
    void
                        initializeRuns(
        Randomness rand);
    // ________________________________________________________________________
    // nat peekNumTax(std::string filePath);
    /**
     *  @brief cleanup, once finished
     */
    void
                        finalizeRuns();
    // ________________________________________________________________________
    /**
     *  @brief starts the MCMC sampling
     */
    void
                        simulate();
    // ________________________________________________________________________
    /**
     *  @brief initializes the config file
     */
    std::tuple<ParameterList, std::vector<std::unique_ptr<AbstractProposal> >,
               std::vector<ProposalSet> >
                        processConfigFile(
        string        configFileName,
        const TreeAln&tralnPtr);
    // ________________________________________________________________________
    void
                        initializeWithParamInitValues(
        TreeAln&            tree,
        const ParameterList&params,
        bool                hasBl) const;
    // ________________________________________________________________________
    void
                        makeTreeUltrametric(
        TreeAln&                        traln,
        std::vector<AbstractParameter*> divTimes,
        std::vector<AbstractParameter*>&divRates) const;
    // ________________________________________________________________________
    /**
     *  @brief EXPERIMENTAL
     */
    void
                        branchLengthsIntegration(
        Randomness&rand);
    // ________________________________________________________________________
    /**
     *  @brief print information about the alignment
     */
    void
                        printAlignmentInfo(
        const TreeAln&traln);
    // ________________________________________________________________________
    /**
     *  @brief gets the runs
     */
    const std::vector<CoupledChains>&
                        getRuns() const
    {return _runs; }
    // ________________________________________________________________________
    virtual void
                        deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void
                        serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    void
                        setCommandLine(
        CommandLine cl){_cl = cl; }
    // ________________________________________________________________________
    void
                        setParallelSetup(
        ParallelSetup* pl){_plPtr = pl; }

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void
                        printParameters(
        const TreeAln&      traln,
        const ParameterList&params) const;
    // ________________________________________________________________________
    void
                        printProposals(
        std::vector<unique_ptr<AbstractProposal> >&proposals,
        std::vector<ProposalSet>&                  proposalSets,
        ParameterList&                             params) const;
    // ________________________________________________________________________
    void
                        printInitializedFiles()
    const;
    // ________________________________________________________________________
    std::vector<std::string>
                        getStartingTreeStrings();
    // ________________________________________________________________________
    void
                        informPrint();
    // ________________________________________________________________________
    void
                        printInitialState();
    // ________________________________________________________________________
    LikelihoodEvaluator
                        createEvaluatorPrototype(
        const TreeAln&initTree,
        bool          useSEV);
    // ________________________________________________________________________
    void
                        writeCheckpointMaster();
    // ________________________________________________________________________
    void
                        initializeFromCheckpoint();
    // ________________________________________________________________________
    std::pair<double,
              double>
                        convergenceDiagnostic(
        nat&begin,
        nat&end);
    // ________________________________________________________________________
    std::vector<bool>
                        initTrees(
        std::vector<TreeAln>&                 trees,
        randCtr_t                             seed,
        std::vector<std::string>              startingTreeStrings,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    bool
                        initializeTree(
        TreeAln&                              traln,
        std::string                           startingTree,
        Randomness&                           treeRandomness,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    void
                        printDuringRun(
        uint64_t gen);
    // ________________________________________________________________________
    std::string
                        getOrCreateBinaryFile()
    const;
    // ________________________________________________________________________
    void
                        catchRunErrors()
    const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<CoupledChains> _runs;
    ParallelSetup*             _plPtr; // non-owning!
    BlockRunParameters         _runParams;
    CommandLine                _cl;
    DiagnosticsFile            _diagFile;

    TimeTracker                _timeTracker;
    TimeTracker                _printTime;
};


#endif
