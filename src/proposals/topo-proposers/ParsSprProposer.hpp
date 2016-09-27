#ifndef _PARS_SPR_DETERMINER_HPP
#define _PARS_SPR_DETERMINER_HPP

#include "TopoMoveProposer.hpp"
#include "SprMove.hpp"
#include "ParsimonyEvaluator.hpp"
#include "Communicator.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             PARS SPR PROPOSER                             //
///////////////////////////////////////////////////////////////////////////////
class ParsSprProposer : public TopoMoveProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using weightMap =  std::unordered_map<BranchPlain, double>;
    using branch2parsScore = std::unordered_map<
                BranchPlain, std::array<parsimonyNumber, 3> >;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
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
    virtual auto                                getMove() const
        ->std::unique_ptr<TopoMove>;
    // ________________________________________________________________________
    virtual void                                printParams(
        std::ostream&out)  const
    {out << ";radius=" << _depth << ",warp=" << _parsWarp;}
    // ________________________________________________________________________
    ParsSprProposer(
        double       parsWarp,
        int          depth,
        Communicator&comm);
    // ________________________________________________________________________
    virtual ~ParsSprProposer(){}
    // ________________________________________________________________________
    ParsSprProposer(
        const ParsSprProposer& rhs)
        : TopoMoveProposer(rhs)
        , _parsWarp{rhs._parsWarp}
        , _depth{rhs._depth}
        , _pEval{rhs._pEval}
        , _comm{rhs._comm}
        , _computedInsertions{rhs._computedInsertions}
        , _move{rhs._move}
    {}
    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    branch2parsScore                            determineScoresOfInsertions(
        TreeAln&    traln,
        BranchPlain prunedTree);
    // ________________________________________________________________________
    void                                        testInsertParsimony(
        TreeAln&         traln,
        nodeptr          insertPos,
        nodeptr          prunedTree,
        branch2parsScore&result,
        int              curDepth);
    // ________________________________________________________________________
    auto                                        getWeights(
        const TreeAln&   traln,
        branch2parsScore insertions)
        ->weightMap;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    double _parsWarp;
    int _depth;
    ParsimonyEvaluator _pEval;
    std::reference_wrapper<Communicator>_comm;
    branch2parsScore _computedInsertions;
    SprMove _move;

    static std::array<double, 2>factors;
    static double weightEps; // mostly a numerical addition
};


#endif
