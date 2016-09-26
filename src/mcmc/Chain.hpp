#ifndef _CHAIN_H
#define  _CHAIN_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include "ProposalSet.hpp"
#include "PriorBelief.hpp"
#include "LikelihoodEvaluator.hpp"
#include "AbstractProposal.hpp"

#include "TopologyFile.hpp"
#include "ParameterFile.hpp"
#include "Serializable.hpp"

#include "CommFlag.hpp"

class TreeAln;
class AbstractProposal;


///////////////////////////////////////////////////////////////////////////////
//                                   CHAIN                                   //
///////////////////////////////////////////////////////////////////////////////
class Chain : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using Param2ContentMap =
            std::unordered_map<AbstractParameter*, ParameterContent>;


    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    Chain(
        randKey_t                                             seed,
        const TreeAln&                                        _traln,
        ParameterList                                         params,
        const std::vector<std::unique_ptr<AbstractProposal> >&_proposals,
        std::vector<ProposalSet>                              proposalSets,
        LikelihoodEvaluator                                   eval,
        bool                                                  isDryRun);
    // ________________________________________________________________________
    Chain(
        const Chain& rhs);
    // ________________________________________________________________________
    Chain(
        Chain&& rhs)  = default;
    // ________________________________________________________________________
    Chain&
                        operator=(
        Chain rhs);

    // void setParamsAndProposals( ,);

    /**
     *  @brief apply saved parameter contents to the tree structure
     *  @param eval indicates whether an evaluation should be performed after
     * resuming
     *  @param checkLnl a hack: disable check for the exact same likelihood.
     * Reason for this is resuming a run from a checkpoint with ExaML. It is
     * just extremely hard to get the exact same likelihood
     */
    // ________________________________________________________________________
    void
                        resume();
    /**
     * @brief saves the all parameters that are integrated over,
     * s.t. the tree can be used by another chain
     * @param paramsOnly indicates whether the likelihood and prior density
     * should be saved as well
     */
    // ________________________________________________________________________
    void
                        suspend();
    /**
     *  @brief proceed by one generation
     */
    // ________________________________________________________________________
    void
                        step();
    /**
     *  @brief gets the proposals of this chain
     */
    // ________________________________________________________________________
    const std::vector<AbstractProposal*>
                        getProposalView()
    const;
    /**
     *  @brief add a representation of the chain to the stream
     */
    // ________________________________________________________________________
    std::ostream&
                        addChainInfo(
        std::ostream&out) const;
    /**
     *  @brief extracts the variables of a chain into a sorted array
     */
    // ________________________________________________________________________
    const ParameterList&
                                            getParameterList()
    const
    {return _params; }
    // ________________________________________________________________________
    /**
     *  @brief take a sample from the chain
     */
    void
                                            sample(
        std::unordered_map<nat, TopologyFile>&tFile,
        ParameterFile&pFile) const;
    // ________________________________________________________________________
    /**
     *  @brief deserialize the input string based on the flags
     */
    void
                                            deserializeConditionally(
        std::istream& in,
        CommFlag      commFlags);
    // ________________________________________________________________________
    /**
     *  @brief serializes the chain into a string based on the flags
     */
    void
                                            serializeConditionally(
        std::ostream& out,
        CommFlag      commFlags) const;
    // ________________________________________________________________________
    void
                                            reinitPrior()
    {_prior.initialize(_traln, _params); }
    // ________________________________________________________________________
    // getters and setters
    log_double
                                            getBestState()
    const
    {return _bestState; }
    // ________________________________________________________________________
    LikelihoodEvaluator&
                                            getEvaluator()
    {return _evaluator; }
    // ________________________________________________________________________
    auto
                                            getTralnHandle()
    const
        ->const TreeAln &
    {return _traln; }

    // ________________________________________________________________________
    auto
                                            getTralnHandle()
        ->TreeAln
    & {return _traln; }
    // ________________________________________________________________________
    Randomness&
                                            getChainRand()
    {return _chainRand; }
    // ________________________________________________________________________
    double
                                            getChainHeat()
    const;
    // ________________________________________________________________________
    void
                                            setDeltaT(
        double dt){_deltaT = dt; }
    // ________________________________________________________________________
    int
                                            getCouplingId()
    const
    {return _couplingId; }
    // ________________________________________________________________________
    void
                                            setCouplingId(
        int id){_couplingId = id; }
    // ________________________________________________________________________
    void
                                            setTuneFreuqency(
        nat _tuneFreq){_tuneFrequency = _tuneFreq; }
    // ________________________________________________________________________
    void
                                            setHeatIncrement(
        nat cplId){_couplingId = cplId; }
    // ________________________________________________________________________
    void
                                            setRunId(
        nat id){_runid = id; }
    // ________________________________________________________________________
    double
                                            getDeltaT()
    {return _deltaT; }
    // ________________________________________________________________________
    uint64_t
                                            getGeneration()
    const
    {return _currentGeneration; }
    // ________________________________________________________________________
    log_double
                                            getLikelihood()
    const
    {return _traln.getLikelihood(); }
    // ________________________________________________________________________
    void
                                            setLikelihood(
        log_double lnl){_traln.setLikelihood(lnl); }
    // ________________________________________________________________________
    void
                                            setLnPr(
        log_double lnPr){_lnPr = lnPr;  }
    // ________________________________________________________________________
    log_double
                                            getLnPr() const
    {return _lnPr; }
    // ________________________________________________________________________
    auto
                                            getPrior()    const
        ->const PriorBelief &
    {return _prior;    }
    // ________________________________________________________________________
    void
                                            updateProposalWeights();
    // ________________________________________________________________________
    const std::vector<ProposalSet>&
                                            getProposalSets()
    const
    {return _proposalSets; }

    // ________________________________________________________________________
    /**
     * @brief the chains determines the virtual root needed by the next
     * generation
     */
    BranchPlain
                                            peekNextVirtualRoot(
        TreeAln&   traln,
        Randomness rand);
    // ________________________________________________________________________
    friend void
                                            swap(
        Chain&lhs,
        Chain&rhs);
    // ________________________________________________________________________
    friend void
                                            swapHeatAndProposals(
        Chain& chainA,
        Chain& chainB);

    // ________________________________________________________________________
    // STUBS
    // currently only here to implement the interface. maybe remove at
    // some point alltogether
    virtual void
                                            deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void
                                            serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    friend std::ostream&
                                            operator
    <<(
        std::ostream& out,
        const Chain&  rhs);
    // ________________________________________________________________________
    void
                                            applyParameterContents(
        Param2ContentMap param2content);
    // ________________________________________________________________________
    // HACK
    void
                                            resetParamPtr();

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    auto
                                            drawProposalFunction(
        Randomness&rand)
        ->AbstractProposal &;
    // ________________________________________________________________________
    auto
                                            drawProposalSet(
        Randomness&rand)
        ->ProposalSet &;
    // ________________________________________________________________________
    void
                                            printArrayStart();
    // ________________________________________________________________________
    void
                                            stepSingleProposal();
    // ________________________________________________________________________
    void
                                            stepSetProposal();

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    TreeAln                                         _traln;

    // this is the global heat parameter that defines the heat
    // increments
    double                                          _deltaT;
    int                                             _runid;
    // TODO should be have per-proposal tuning?
    int                                             _tuneFrequency;
    log_double                                      _hastings;         // logged!
    uint64_t                                        _currentGeneration;
    /// indicates how hot the chain is (i = 0 => cold chain), may change!
    nat                                             _couplingId;         // CHECKPOINTED
    std::vector<std::unique_ptr<AbstractProposal> > _proposals;
    std::vector<ProposalSet>                        _proposalSets;
    Randomness                                      _chainRand;
    double                                          _relWeightSumSingle;
    double                                          _relWeightSumSets;
    PriorBelief                                     _prior;
    log_double                                      _bestState;
    LikelihoodEvaluator                             _evaluator;

    // suspending and resuming the chain
    log_double                                      _lnPr;

    ParameterList                                   _params;
};


#endif
