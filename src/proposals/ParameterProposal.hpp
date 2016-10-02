#ifndef _PARAMETER_PROPOSALS
#define _PARAMETER_PROPOSALS

#include "TreeAln.hpp"
#include "AbstractProposal.hpp"
#include "AbstractProposer.hpp"
#include "ParameterProposal.hpp"

#include <memory>

///////////////////////////////////////////////////////////////////////////////
//                             PARAMETER PROPOSAL                            //
///////////////////////////////////////////////////////////////////////////////
class ParameterProposal : public AbstractProposal
{
public:
    // ________________________________________________________________________
    ParameterProposal(
        Category               cat,
        std::string            name,
        bool                   modifiesBL,
        AbstractProposer::UPtr proposer,
        double                 parameter,
        double                 weight,
        double                 minTuning,
        double                 maxTuning);
    // ________________________________________________________________________
    ParameterProposal(
        const ParameterProposal& prop);
    // ________________________________________________________________________
    virtual ~ParameterProposal(){}
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
    virtual BranchPlain                               determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const {return BranchPlain();}
    // ________________________________________________________________________
    virtual AbstractProposal*                         clone()
    const {return new ParameterProposal(*this);}
    // ________________________________________________________________________
    virtual std::vector<nat>                          getInvalidatedNodes(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual std::pair<BranchPlain,
                      BranchPlain>                    prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand)
    {return std::make_pair(BranchPlain(0, 0), BranchPlain(0, 0));}
    // ________________________________________________________________________
    virtual void                                      readFromCheckpointCore(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                      writeToCheckpointCore(
        std::ostream&out) const;
private:
    bool modifiesBL;
    double parameter;
    AbstractProposer::UPtr proposer;

    ParameterContent _savedContent;
    ParameterContent _savedBinaryContent;
    std::vector<double>_oldFCs;
};


#endif
