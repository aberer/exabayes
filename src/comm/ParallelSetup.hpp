#ifndef _PARALLEL_SETUP_PARA_HPP
#define _PARALLEL_SETUP_PARA_HPP

#include "CommandLine.hpp"
#include "IncompleteMesh.hpp"
#include "FlagType.hpp"
#include "Communicator.hpp"

#include "ThreadResource.hpp"

#include "CommFlag.hpp"

#include <unordered_map>

class CoupledChains;

// remember: the implementation of this class depends on whether we
// have MPI. Do NOT implement anything here inline.


///////////////////////////////////////////////////////////////////////////////
//                               PARALLEL SETUP                              //
///////////////////////////////////////////////////////////////////////////////
class ParallelSetup
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ParallelSetup(
        CommandLine& cl);
    // ________________________________________________________________________
    ParallelSetup(
        ParallelSetup&&rhs) = default;
    // ________________________________________________________________________
    ParallelSetup(
        const ParallelSetup& rhs) = delete;
    // ________________________________________________________________________
    ParallelSetup&                         operator=(
        ParallelSetup rhs);
    // ________________________________________________________________________
    friend void                            swap(
        ParallelSetup&lhs,
        ParallelSetup&rhs);
    // ________________________________________________________________________
    void                                   warn();
    // ________________________________________________________________________
    /**
     *  @brief initializes the examl environment
     */
    void                                   initialize();
    // ________________________________________________________________________
    static void                            initializeRemoteComm(
        int   argc,
        char**argv);
    // ________________________________________________________________________
    /**
     *  @brief finalize the parallel environment
     */
    static void                            finalize();
    // ________________________________________________________________________
    /**
     *  @brief gets the number of runs executed in parallel
     */
    size_t                                 getRunsParallel() const;
    // ________________________________________________________________________
    /**
     *  @brief gets the number of (coupled) chains executed in paralel (in
     * addition to run-level parallelism)
     */
    // ________________________________________________________________________
    size_t                                 getChainsParallel() const;
    /**
     * @brief indicates whether the process should conduct output operatons
     * (for his run)
     */
    // ________________________________________________________________________
    bool                                   isRunLeader() const;
    // ________________________________________________________________________
    bool                                   isRunLeader(
        nat gRank) const;
    // ________________________________________________________________________
    /**
     * @brief indicates whether a process is the leader in a batch of chains
     */
    bool                                   isChainLeader() const;
    // ________________________________________________________________________
    /**
     *  @brief indicates whether the process is the process to conduct global
     * logging output
     */
    bool                                   isGlobalMaster() const;
    // ________________________________________________________________________
    /**
     *  @brief sends a serialized chain representation to the other involved
     * process
     *  @param run relevant run
     *  @param otherChainIndex the index of the chain that is involved in the
     * swapping attempt
     *  @return serialized representation  of the chain we swap with
     */
    // ________________________________________________________________________
    std::string                            sendRecvChain(
        const CoupledChains& run,
        nat                  myIndex,
        nat                  otherChainIndex,
        std::string          myChainSer,
        CommFlag             flags);
    // ________________________________________________________________________
    nat                                    chainIdToLeaderRank(
        nat runId,
        nat chainId) const;
    // ________________________________________________________________________
    /**
     *  @brief synchronize all core chain information at the master node //
     *  @notice this only concerns the global masterprocess in a
     *  parallel environtment.
     */
    void                                   synchronizeChainsAtMaster(
        std::vector<CoupledChains>& runs,
        CommFlag                    commFlags);
    // ________________________________________________________________________
    /**
     *  @brief a printing function that allows to print properly
     */
    static void                            blockingPrint(
        Communicator&comm,
        std::string  ss);
    // ________________________________________________________________________
    std::string                            printLoadBalance(
        const TreeAln& traln,
        nat            numRun,
        nat            numCoupled);
    // ________________________________________________________________________
    /**
     *  @brief indicates whether the chain belongs to a process
     */
    bool                                   isMyChain(
        nat runid,
        nat chainIndex) const;
    // ________________________________________________________________________
    bool                                   isMyChain(
        nat gRank,
        nat runid,
        nat chainIndex) const;
    // ________________________________________________________________________
    /**
     *  @brief indicates whether the run belongs to a process
     */
    bool                                   isMyRun(
        nat runid) const;
    // ________________________________________________________________________
    bool                                   isMyRun(
        nat gRank,
        nat runid) const;
    // ________________________________________________________________________
    Communicator&                          getGlobalComm();
    // ________________________________________________________________________
    Communicator&                          getRunComm();
    // ________________________________________________________________________
    Communicator&                          getChainComm();
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&        out,
        const ParallelSetup& rhs);
    // ________________________________________________________________________
    IncompleteMesh                         getMesh() const;
    // ________________________________________________________________________
    int                                    getRankInRun(
        std::array<nat, 3>coords) const;
    // ________________________________________________________________________
    std::array<nat, 3>                     getMyCoordinates() const;
    // ________________________________________________________________________
    static void                            abort(
        int  code,
        bool waitForAll);
    // ________________________________________________________________________
    void                                   releaseThreads();
    // ________________________________________________________________________
    /**
     * @brief determines, whether the swap can be executed within a
     * local communicator only
     */
    bool                                   swapIsLocal(
        nat chainIdA,
        nat chainIdB,
        nat runId) const;
    // ________________________________________________________________________
    uint64_t                               getMaxTag();
    // ________________________________________________________________________
    void                                   pinThreads();

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    auto                                   serializeAllChains(
        std::vector<CoupledChains>&runs,
        CommFlag                   commFlags) const
        ->std::tuple<nat, std::vector<char> >;
    // ________________________________________________________________________
    auto                                   serializeSwapMatrices(
        std::vector<CoupledChains>& runs,
        CommFlag&                   commFlags) const
        ->std::tuple<nat, std::vector<char> >;
    // ________________________________________________________________________
    int                                    getColorInChain(
        std::array<nat, 3>coords) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    ThreadResource _threadResource;

    Communicator   _globalComm;
    Communicator   _runComm;
    Communicator   _chainComm;

    IncompleteMesh _mesh;
};


#endif
