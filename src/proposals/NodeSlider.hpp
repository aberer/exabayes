#ifndef _NODESLIDER_H
#define _NODESLIDER_H

#include "AbstractProposal.hpp"
#include "Path.hpp"

///////////////////////////////////////////////////////////////////////////////
//                                NODE SLIDER                                //
///////////////////////////////////////////////////////////////////////////////
class NodeSlider : public AbstractProposal
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    NodeSlider(
        double multiplier);
    // ________________________________________________________________________
    virtual ~NodeSlider(){}
    // ________________________________________________________________________
    virtual auto                                determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const
        ->BranchPlain;
    // ________________________________________________________________________
    virtual void                                applyToState(
        TreeAln&             traln,
        PriorBelief&         prior,
        log_double&          hastings,
        Randomness&          rand,
        LikelihoodEvaluator& eval);
    // ________________________________________________________________________
    virtual void                                evaluateProposal(
        LikelihoodEvaluator&evaluator,
        TreeAln&            traln,
        const BranchPlain&  branchSuggestion);
    // ________________________________________________________________________
    virtual void                                resetState(
        TreeAln&traln);
    // ________________________________________________________________________
    BranchPlain                                 proposeBranch(
        const TreeAln&traln,
        Randomness&   rand) const;
    // ________________________________________________________________________
    BranchPlain                                 proposeOtherBranch(
        const BranchPlain&firstBranch,
        const TreeAln&    traln,
        Randomness&       rand) const;
    // ________________________________________________________________________
    virtual void                                prepareForSetEvaluation(
        TreeAln&             traln,
        LikelihoodEvaluator& eval) const;
    // ________________________________________________________________________
    auto                                        prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand)
        ->std::pair<BranchPlain, BranchPlain>;
    // ________________________________________________________________________
    virtual void                                autotune()
    {}                          // disabled
    // ________________________________________________________________________
    virtual AbstractProposal*                   clone() const;
    // ________________________________________________________________________
    virtual auto                                getInvalidatedNodes(
        const TreeAln&traln) const
        ->std::vector<nat>;
    // ________________________________________________________________________
    virtual void                                readFromCheckpointCore(
        std::istream&in){}                                    // disabled
    // ________________________________________________________________________
    virtual void                                writeToCheckpointCore(
        std::ostream&out) const {}                                // disabled

    ///////////////////////////////////////////////////////////////////////////
    //                            PROTECTED DATA                             //
    ///////////////////////////////////////////////////////////////////////////
protected:
    double multiplier;
    BranchLength oneBranch;
    BranchLength otherBranch;
};


#endif
