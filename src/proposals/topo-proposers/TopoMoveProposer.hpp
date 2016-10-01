#ifndef _ABSTRACT_TOPO_MOVE_DRAWER_HPP
#define _ABSTRACT_TOPO_MOVE_DRAWER_HPP

#include "TreeAln.hpp"
#include "log_double.hpp"
#include "TopoMove.hpp"
#include "LikelihoodEvaluator.hpp"
#include "OptimizedParameter.hpp"

#include <memory>


class Communicator;

///////////////////////////////////////////////////////////////////////////////
//                             TOPO MOVE PROPOSER                            //
///////////////////////////////////////////////////////////////////////////////
class TopoMoveProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using UPtr = std::unique_ptr<TopoMoveProposer>;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    /**
     *  determines the move and a probability of drawing this move
     */
    virtual void                                         determineMove(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        Randomness&                           rand,
        BranchPlain                           primeBranch,
        const std::vector<AbstractParameter*>&params)  = 0;
    // ________________________________________________________________________
    /**
     *  only determines the proposal density of the backward move (NOT
     * including any further transformations like branch lengths)
     */
    virtual void                                         determineBackProb(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        const std::vector<AbstractParameter*>&params) = 0;
    // ________________________________________________________________________
    /**
     * determines the first branch (e.g., subtree that is pruned)
     */
    virtual BranchPlain                                  determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const = 0;
    // ________________________________________________________________________
    virtual TopoMoveProposer*                            clone() const = 0;
    // ________________________________________________________________________
    virtual std::unique_ptr<TopoMove>                    getMove() const  = 0;
    // ________________________________________________________________________
    virtual void                                         printParams(
        std::ostream&out)  const {}
    // ________________________________________________________________________
    TopoMoveProposer()
        : _forwProb{log_double::fromAbs(1.)}
        , _backProb{log_double::fromAbs(1.)}
        , _withOptimizedBranches{false}
        , _optBranches{}
    {}
    // ________________________________________________________________________
    virtual ~TopoMoveProposer(){}
    // ________________________________________________________________________
    // ________________________________________________________________________
    log_double                                           getProposalDensity()
    const {return _backProb / _forwProb;}
    // ________________________________________________________________________
    bool                                                 hasOptimizedBranches()
    const {return _withOptimizedBranches;}
    // ________________________________________________________________________
    std::vector<OptimizedParameter>                      getOptimizedBranches()
    const {return _optBranches;}

    ///////////////////////////////////////////////////////////////////////////
    //                            PROTECTED DATA                             //
    ///////////////////////////////////////////////////////////////////////////
protected:
    log_double _forwProb;
    log_double _backProb;
    bool _withOptimizedBranches;
    std::vector<OptimizedParameter>_optBranches;
};


#endif
