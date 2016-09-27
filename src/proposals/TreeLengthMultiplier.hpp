#include "AbstractProposal.hpp"


///////////////////////////////////////////////////////////////////////////////
//                           TREE LENGTH MULTIPLIER                          //
///////////////////////////////////////////////////////////////////////////////
class TreeLengthMultiplier : public AbstractProposal
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    TreeLengthMultiplier(
        double _multiplier);
    // ________________________________________________________________________
    virtual ~TreeLengthMultiplier(){}
    // ________________________________________________________________________
    virtual auto                    determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const
        ->BranchPlain
    {return BranchPlain();}
    // ________________________________________________________________________
    virtual void                    applyToState(
        TreeAln&             traln,
        PriorBelief&         prior,
        log_double&          hastings,
        Randomness&          rand,
        LikelihoodEvaluator& eval);
    // ________________________________________________________________________
    virtual void                    evaluateProposal(
        LikelihoodEvaluator&evaluator,
        TreeAln&            traln,
        const BranchPlain&  branchSuggestion);
    // ________________________________________________________________________
    virtual void                    resetState(
        TreeAln&traln);
    // ________________________________________________________________________
    virtual void                    autotune();
    // ________________________________________________________________________
    virtual auto                    clone()
    const
        ->AbstractProposal *;
    // ________________________________________________________________________
    virtual void                    readFromCheckpointCore(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                    writeToCheckpointCore(
        std::ostream&out) const;
    // ________________________________________________________________________
    virtual auto                    prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand)
        ->std::pair<BranchPlain, BranchPlain>
    {return std::make_pair(BranchPlain(0, 0), BranchPlain(0, 0));}
    // ________________________________________________________________________
    virtual auto                    getInvalidatedNodes(
        const TreeAln& traln) const
        ->std::vector<nat>;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                            multiplyBranchLengthsRecursively(
        TreeAln& traln,
        nodeptr  p,
        double   multiHere);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double multiplier;      // the tuning variable
    double initTreeLength;
    std::vector<BranchLength>storedBranches;
};

