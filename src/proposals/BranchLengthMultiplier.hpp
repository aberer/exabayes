#ifndef _BRANCHLENGTHSMULTIPLIER_H
#define _BRANCHLENGTHSMULTIPLIER_H

#include "AbstractProposal.hpp"

#include <limits>


///////////////////////////////////////////////////////////////////////////////
//                          BRANCH LENGTH MULTIPLIER                         //
///////////////////////////////////////////////////////////////////////////////
class BranchLengthMultiplier : public AbstractProposal
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BranchLengthMultiplier(
        double multiplier);
    // ________________________________________________________________________
    virtual BranchPlain                               determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const;
    // ________________________________________________________________________
    virtual void                                      applyToState(
        TreeAln&             traln,
        PriorBelief&         prior,
        log_double&          hastings,
        Randomness&          rand,
        LikelihoodEvaluator& eval);
    // ________________________________________________________________________
    virtual void                                      evaluateProposal(
        LikelihoodEvaluator&evaluator,
        TreeAln&            traln,
        const BranchPlain&  branchSuggestion);
    // ________________________________________________________________________
    virtual void                                      resetState(
        TreeAln&traln);
    // ________________________________________________________________________
    virtual void                                      autotune();
    // ________________________________________________________________________
    virtual std::vector<nat>                          getInvalidatedNodes(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual AbstractProposal*                         clone() const;
    // ________________________________________________________________________
    virtual std::pair<BranchPlain,
                      BranchPlain>                    prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand);
    // ________________________________________________________________________
    virtual BranchPlain                               proposeBranch(
        const TreeAln&traln,
        Randomness&   rand) const;
    // ________________________________________________________________________
    virtual void                                      readFromCheckpointCore(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                      writeToCheckpointCore(
        std::ostream&out) const;


    ///////////////////////////////////////////////////////////////////////////
    //                             PROTECTED DATA                            //
    ///////////////////////////////////////////////////////////////////////////
protected:
    double       _multiplier;
    BranchLength _savedBranch;
};


#endif
