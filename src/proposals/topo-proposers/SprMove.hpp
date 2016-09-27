/**
 *  @brief a helper object making the code concerning a SPR move
 *  re-usable.
 */


#ifndef _SPR_MOVE_H
#define _SPR_MOVE_H

#include "Path.hpp"

#include "TreeTraverser.hpp"

#include "TopoMove.hpp"

#include <unordered_set>

class PriorBelief;
class LikelihoodEvaluator;

///////////////////////////////////////////////////////////////////////////////
//                                  SPR MOVE                                 //
///////////////////////////////////////////////////////////////////////////////
class SprMove : public TopoMove
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void                                         apply(
        TreeAln&                              traln,
        const std::vector<AbstractParameter*>&blParams) const;
    // ________________________________________________________________________
    virtual BranchPlain                                  getEvalBranch(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual auto                                         getInverse() const
        ->std::unique_ptr<TopoMove>;
    // ________________________________________________________________________
    virtual TopoMove*                                    clone() const;
    // ________________________________________________________________________
    virtual void                                         print(
        std::ostream& out) const;
    // ________________________________________________________________________
    virtual std::vector<BranchPlain>                     getBranchesToPropose(
        const TreeAln& traln,
        MoveOptMode    mode);
    // ________________________________________________________________________
    virtual int                                          getNniDistance() const
    {return std::max(int(_path.size()) - 2, 0);}
    // ________________________________________________________________________
    virtual std::vector<nat>                             getDirtyNodes() const;

    virtual void                                         invalidateArrays(
        LikelihoodEvaluator& eval,
        TreeAln&             traln,
        MoveOptMode          mode)  const;
    // ________________________________________________________________________
    explicit SprMove()
        : _path{}
    {}
    // ________________________________________________________________________
    virtual ~SprMove(){}
    // ________________________________________________________________________
    SprMove(
        const TreeAln&    traln,
        const BranchPlain&prunedTree,
        const BranchPlain&insertBranch);
    // ________________________________________________________________________
    auto                                                 getDirtyWrtPrevious(
        const SprMove&rhs) const
        ->std::vector<nat>;
    // ________________________________________________________________________
    BranchPlain                                          getMovingSubtree(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    auto                                                 getNodeList() const
        ->std::vector<nat>;
    // ________________________________________________________________________
    SprMove                                              getInverseMove() const;
    // ________________________________________________________________________
    BranchPlain                                          getInsertBranch()
    const
    {return _path.at(int(_path.size()) - 1);}
    // ________________________________________________________________________
    auto                                                 getBranchesByDistance(
        const TreeAln&traln,
        nat           dist) const
        ->std::unordered_set<BranchPlain>;
    // ________________________________________________________________________
    friend std::ostream&                                 operator<<(
        std::ostream&  out,
        const SprMove& rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    Path _path;
};


#endif
