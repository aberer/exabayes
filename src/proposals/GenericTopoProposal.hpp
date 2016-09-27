#ifndef _GENERIC_TOPO_PROPOSAL_HPP
#define _GENERIC_TOPO_PROPOSAL_HPP

#include "AbstractProposal.hpp"
#include "TopoMoveProposer.hpp"
#include "BranchBackup.hpp"

///////////////////////////////////////////////////////////////////////////////
//                           GENERIC TOPO PROPOSAL                           //
///////////////////////////////////////////////////////////////////////////////
class GenericTopoProposal : public AbstractProposal
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void
                                                                                                                  applyToState(
        TreeAln&             traln,
        PriorBelief&         prior,
        log_double&          hastings,
        Randomness&          rand,
        LikelihoodEvaluator& eval);
    // ________________________________________________________________________
    virtual void
                                                                                                                  evaluateProposal(
        LikelihoodEvaluator&evaluator,
        TreeAln&            traln,
        const BranchPlain&  branchSuggestion);
    // ________________________________________________________________________
    virtual void
                                                                                                                  resetState(
        TreeAln&traln);
    virtual std::pair<BranchPlain,
                      BranchPlain>
                                                                                                                  prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand)
    {return std::make_pair(BranchPlain(0, 0), BranchPlain(0, 0));}
    // ________________________________________________________________________
    virtual void
                                                                                                                  autotune()
    {}
    // ________________________________________________________________________
    virtual auto
    clone()
    const
        ->AbstractProposal *;
    // ________________________________________________________________________
    virtual void
                                                                                                                  readFromCheckpointCore(
        std::istream&in)
    {}
    // ________________________________________________________________________
    virtual void
                                                                                                                  writeToCheckpointCore(
        std::ostream&out) const {}                                // disabled
    // ________________________________________________________________________
    virtual void
                                                                                                                  printParams(
        std::ostream&out)  const {_moveProposer->printParams(out);}

    // ________________________________________________________________________
    GenericTopoProposal(
        std::unique_ptr<TopoMoveProposer>moveDet,
        std::string                      name,
        double                           relWeight,
        MoveOptMode                      toOpt);
    // ________________________________________________________________________
    GenericTopoProposal(
        const GenericTopoProposal& rhs);
    // ________________________________________________________________________
    GenericTopoProposal(
        GenericTopoProposal&& rhs) = delete;
    // ________________________________________________________________________
    GenericTopoProposal&
                                                                                                                  operator
    =(
        GenericTopoProposal rhs);
    // ________________________________________________________________________
    friend void
    swap(
        GenericTopoProposal& lhs,
        GenericTopoProposal& rhs);
    // ________________________________________________________________________
    void
                                                                                                                  enableUseMultiplier()
    {_useMultiplier = true;}
    // ________________________________________________________________________
    log_double
                                                                                                                  assessOldBranches(
        TreeAln&                             traln,
        LikelihoodEvaluator&                 eval,
        const std::vector<BranchPlain>&      bs,
        const std::vector<AbstractParameter*>params);
    // ________________________________________________________________________
    BranchPlain
                                                                                                                  determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const;
    // ________________________________________________________________________
    std::vector<nat>
                                                                                                                  getInvalidatedNodes(
        const TreeAln& traln) const;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    log_double
                                                                                                                  proposeAndApplyNewBranches(
        TreeAln&                       traln,
        const std::vector<BranchPlain>&bs,
        PriorBelief&                   prior,
        Randomness&                    rand,
        LikelihoodEvaluator&           eval) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::unique_ptr<TopoMove>_move; // abstract product
    std::unique_ptr<TopoMoveProposer>_moveProposer; // abstract factory
    MoveOptMode _moveOptMode;
    BranchBackup _backup;
    bool _useMultiplier;

    static double multiplier;
};


#endif
