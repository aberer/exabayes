#ifndef LIKELIHOODSPRPROPOSER_H
#define LIKELIHOODSPRPROPOSER_H

#include "TopoMoveProposer.hpp"
#include "SprMove.hpp"
#include "BranchBackup.hpp"

///////////////////////////////////////////////////////////////////////////////
//                           LIKEHOOD SPR PROPOSER                           //
///////////////////////////////////////////////////////////////////////////////
class LikehoodSprProposer : public TopoMoveProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void                              determineMove(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        Randomness&                           rand,
        BranchPlain                           primeBranch,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    virtual void                              determineBackProb(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    virtual auto                              determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const
        ->BranchPlain;
    // ________________________________________________________________________
    virtual auto                              clone() const
        ->TopoMoveProposer *;
    // ________________________________________________________________________
    virtual TopoMove::UPtr                    getMove() const;
    // ________________________________________________________________________
    virtual void                              printParams(
        std::ostream&out)  const
    {out << ";radius=" << _maxStep << ",warp=" << _likeWarp;}

    // ________________________________________________________________________
    LikehoodSprProposer(
        nat         maxStep,
        double      likeWarp,
        MoveOptMode toOpt);
    virtual ~LikehoodSprProposer(){}

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    auto                                      computeLikelihoodsOfInsertions(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        const BranchPlain&                    prunedTree,
        const std::vector<AbstractParameter*>&params)
        ->std::vector<InsertionResult>;
    // ________________________________________________________________________
    auto                                      transformLikelihoods(
        std::vector<InsertionResult>result) const
        ->std::vector<InsertionResult>;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    nat _maxStep;
    double _likeWarp;
    SprMove _move;
    MoveOptMode _moveOptMode;

    static double weightEps;
};


#endif /* LIKELIHOODSPRPROPOSER_H */
