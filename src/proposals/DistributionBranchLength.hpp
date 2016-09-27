#ifndef _DISTRBUTION_BRANCH_LENGTH
#define  _DISTRBUTION_BRANCH_LENGTH

#include "BranchLengthMultiplier.hpp"

#include "DistributionProposer.hpp"
#include "GammaProposer.hpp"

#include "ComplexTuner.hpp"

#include "OptimizedParameter.hpp"

///////////////////////////////////////////////////////////////////////////////
//                         DISTRIBUTION BRANCH LENGTH                        //
///////////////////////////////////////////////////////////////////////////////
template<class C>
class DistributionBranchLength : public BranchLengthMultiplier
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    DistributionBranchLength();
    // ________________________________________________________________________
    virtual void                                      createProposer(
        TreeAln&             traln,
        LikelihoodEvaluator& eval,
        BranchPlain          b);
    // ________________________________________________________________________
    virtual void                                      applyToState(
        TreeAln&             traln,
        PriorBelief&         prior,
        log_double&          hastings,
        Randomness&          rand,
        LikelihoodEvaluator& eval);
    // ________________________________________________________________________
    virtual void                                      autotune();
    // ________________________________________________________________________
    virtual std::pair<BranchPlain,
                      BranchPlain>                    prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand);
    // ________________________________________________________________________
    virtual AbstractProposal*                         clone() const
    {return new DistributionBranchLength(*this);}
    // ________________________________________________________________________
    virtual void                                      extractProposer(
        TreeAln&                  traln,
        const OptimizedParameter& param)
    {
        auto attr = _allParams->at(_primParamIds[0])->getAttributes();
        _proposer = param.getProposerDistribution<C>(traln,
                                                     attr._convTuner.
                                                         getParameter(),
                                                     attr._nonConvTuner.
                                                         getParameter());
    }
    // ________________________________________________________________________
    virtual void                                      accept();
    // ________________________________________________________________________
    virtual void                                      reject();

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    DistributionProposer<C>_proposer;
    bool _converged;
};


// ____________________________________________________________________________
void                                                  calcDiagptable(
    const double z,
    const int    states,
    const int    numberOfCategories,
    const double*rptr,
    const double*EIGN,
    double*      diagptable);

#include "DistributionBranchLength.tpp"

#endif
