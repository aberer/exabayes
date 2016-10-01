#ifndef _EXTENDED_TBR_DETERMINER_HPP
#define _EXTENDED_TBR_DETERMINER_HPP

#include "TbrMove.hpp"
#include "TopoMoveProposer.hpp"

class Path;

///////////////////////////////////////////////////////////////////////////////
//                           EXTENDED TBR PROPOSER                           //
///////////////////////////////////////////////////////////////////////////////
class ExtendedTbrProposer : public TopoMoveProposer
{
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
    TopoMove::UPtr                              getMove() const;
    // ________________________________________________________________________
    virtual void                                printParams(
        std::ostream&out)  const {out << ";stopProb=" << _stopProb;}
    // ________________________________________________________________________
    ExtendedTbrProposer(
        double stopProb)
        : _stopProb{stopProb}
        , _move{}
    {}

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                                        buildPath(
        Path&                          path,
        BranchPlain                    bisectedBranch,
        TreeAln&                       traln,
        Randomness&                    rand,
        std::vector<AbstractParameter*>params) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    double _stopProb;
    TbrMove _move;
};


#endif
