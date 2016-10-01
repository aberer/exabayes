#ifndef _TBR_MOVE_H
#define _TBR_MOVE_H

#include "SprMove.hpp"

#include "TopoMove.hpp"

///////////////////////////////////////////////////////////////////////////////
//                                  TBR MOVE                                 //
///////////////////////////////////////////////////////////////////////////////
class TbrMove : public TopoMove
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using UPtr = std::unique_ptr<TbrMove>;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void                                         apply(
        TreeAln&                              traln,
        const std::vector<AbstractParameter*>&blParams) const;
    // ________________________________________________________________________
    virtual ~TbrMove(){}

    virtual BranchPlain                                  getEvalBranch(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual std::vector<nat>                             getDirtyNodes() const;
    // ________________________________________________________________________
    virtual TopoMove::UPtr                               getInverse() const;
    // ________________________________________________________________________
    virtual TopoMove*                                    clone() const;
    // ________________________________________________________________________
    virtual void                                         print(
        std::ostream& out) const;
    // ________________________________________________________________________
    virtual auto                                         getBranchesToPropose(
        const TreeAln& traln,
        MoveOptMode    mode)
        ->std::vector<BranchPlain>;
    // ________________________________________________________________________
    virtual int                                          getNniDistance() const
    {return _moveA.getNniDistance() + _moveB.getNniDistance();}
    // ________________________________________________________________________
    virtual void                                         invalidateArrays(
        LikelihoodEvaluator& eval,
        TreeAln&             traln,
        MoveOptMode          mode) const;
    // ________________________________________________________________________
    explicit TbrMove()
        : _moveA{}
        , _moveB{}
    {}
    // ________________________________________________________________________
    TbrMove(
        const TreeAln&     traln,
        const BranchPlain& bisected,
        const BranchPlain& insertA,
        const BranchPlain& insertB);
    // ________________________________________________________________________
    TbrMove                                              getInverseMove() const;
    // ________________________________________________________________________
    friend std::ostream&                                 operator<<(
        std::ostream& out,
        const TbrMove&rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    SprMove _moveA;
    SprMove _moveB;
};


#endif
