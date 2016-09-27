#ifndef _STAT_NNI_HPP_NEW
#define _STAT_NNI_HPP_NEW

#include "SprMove.hpp"

#include "TopoMoveProposer.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             STAT NNI PROPOSER                             //
///////////////////////////////////////////////////////////////////////////////
class StatNniProposer : public TopoMoveProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    StatNniProposer()
        : _move{}
    {}
    // ________________________________________________________________________
    virtual void                                determineMove(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        Randomness&                           rand,
        BranchPlain                           primeBranch,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    virtual TopoMoveProposer*                   clone() const;
    // ________________________________________________________________________
    virtual BranchPlain                         determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const;
    // ________________________________________________________________________
    virtual void                                determineBackProb(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    auto                                        getMove() const
        ->std::unique_ptr<TopoMove>;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    SprMove _move;
};


#endif
