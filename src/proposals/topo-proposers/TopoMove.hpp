#ifndef _ABSTRACT_MOVE_HPP
#define _ABSTRACT_MOVE_HPP

#include "TreeAln.hpp"
#include "LikelihoodEvaluator.hpp"

#include "MoveOptMode.hpp"

#include <vector>

class AbstractParameter;

///////////////////////////////////////////////////////////////////////////////
//                                 TOPO MOVE                                 //
///////////////////////////////////////////////////////////////////////////////
class TopoMove
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using UPtr = std::unique_ptr<TopoMove>;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual ~TopoMove()
    {}
    // ________________________________________________________________________
    virtual void                                         apply(
        TreeAln&                              traln,
        const std::vector<AbstractParameter*>&blParams) const = 0;
    // ________________________________________________________________________
    virtual BranchPlain                                  getEvalBranch(
        const TreeAln&traln) const = 0;
    // ________________________________________________________________________
    virtual std::vector<nat>                             getDirtyNodes() const
        =
            0;
    // ________________________________________________________________________
    virtual auto                                         getInverse() const
        ->UPtr = 0;
    // ________________________________________________________________________
    virtual TopoMove*                                    clone() const = 0;
    // ________________________________________________________________________
    virtual void                                         print(
        std::ostream& out) const = 0;
    // ________________________________________________________________________
    virtual int                                          getNniDistance() const
        =
            0;
    // ________________________________________________________________________
    virtual std::vector<BranchPlain>                     getBranchesToPropose(
        const TreeAln& traln,
        MoveOptMode    mode) = 0;
    // ________________________________________________________________________
    virtual void                                         invalidateArrays(
        LikelihoodEvaluator& eval,
        TreeAln&             traln,
        MoveOptMode          mode) const = 0;
    // ________________________________________________________________________
    friend std::ostream&                                 operator<<(
        std::ostream&   out,
        const TopoMove& elem)
    {
        elem.print(out);
        return out;
    }
};


#endif
