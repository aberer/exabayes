#include "TopoMoveProposer.hpp"

#include "SprMove.hpp"

///////////////////////////////////////////////////////////////////////////////
//                           EXTENDED SPR PROPOSER                           //
///////////////////////////////////////////////////////////////////////////////
class ExtendedSprProposer : public TopoMoveProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void                                         determineMove(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        Randomness&                           rand,
        BranchPlain                           primeBranch,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    virtual TopoMoveProposer*                            clone() const;
    // ________________________________________________________________________
    virtual BranchPlain                                  determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const;
    // ________________________________________________________________________
    virtual void                                         determineBackProb(
        TreeAln&                              traln,
        LikelihoodEvaluator&                  eval,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    virtual std::unique_ptr<TopoMove>                    getMove() const;
    // ________________________________________________________________________
    virtual void                                         printParams(
        std::ostream&out)  const {out << ";stopProb=" << _stopProb;}
    // ________________________________________________________________________
    ExtendedSprProposer(
        double stopProb);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double _stopProb;
    SprMove _move;
};

