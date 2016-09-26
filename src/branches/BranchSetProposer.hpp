#ifndef BRANCH_SET_PROPOSER_HPP
#define BRANCH_SET_PROPOSER_HPP

#include <unordered_map>

#include "OptimizedParameter.hpp"

class LikelihoodEvaluator;


///////////////////////////////////////////////////////////////////////////////
//                            BRANCH SET PROPOSER                            //
///////////////////////////////////////////////////////////////////////////////
class BranchSetProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
    using ResultType =
            std::unordered_map<BranchPlain,
                               std::vector<OptimizedParameter> >;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ~BranchSetProposer(){}
    // ________________________________________________________________________
    BranchSetProposer(
        TreeAln&                       traln,
        std::vector<BranchPlain>       branches,
        std::vector<AbstractParameter*>params);
    // ________________________________________________________________________
    void                          findJointOptimum(
        LikelihoodEvaluator& eval,
        int                  maxIter,
        bool                 computeLikelihood);
    // ________________________________________________________________________

    auto                          getResult() const
        ->ResultType
    {return _result; }
    // ________________________________________________________________________
    log_double                    getOptimalLikelihood()
    const
    {
        return _likelihood;
    }
    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                          reorderToConnectedComponent();

    ///////////////////////////////////////////////////////////////////////////
    //                             PRIVATE DATA                              //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::reference_wrapper<TreeAln>                                   _traln;
    std::vector<BranchPlain>                                          _branches;
    const std::vector<AbstractParameter*>                             _params;
    std::unordered_map<BranchPlain, std::vector<OptimizedParameter> > _result;
    log_double _likelihood;
};


#endif
