#ifndef BRANCHBASE_H
#define BRANCHBASE_H

#include "common.h"

#include "Serializable.hpp"

#include <cassert>

#include <fstream>

///////////////////////////////////////////////////////////////////////////////
//                                BRANCH PLAIN                               //
///////////////////////////////////////////////////////////////////////////////
class BranchPlain : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    explicit BranchPlain(
        nat a = 0,
        nat b = 0)
        : _primNode{a}
        , _secNode{b}
    {}

    // implements Serializable
    // ________________________________________________________________________
    virtual void                           deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                           serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    virtual ~BranchPlain(){}
    // ________________________________________________________________________
    BranchPlain                            getInverted() const
    {return BranchPlain(_secNode, _primNode); }
    // ________________________________________________________________________
    bool                                   operator==(
        const BranchPlain&other) const
    {return other._primNode == _primNode && other._secNode == _secNode; }
    // ________________________________________________________________________
    bool                                   operator!=(
        const BranchPlain&other) const
    {return not (*this == other); }
    // ________________________________________________________________________
    bool                                   hasNode(
        nat aNode) const
    {return _primNode == aNode  || _secNode == aNode; }
    // ________________________________________________________________________
    nat                                    getPrimNode() const
    {return _primNode; }
    // ________________________________________________________________________
    nat                                    getSecNode() const
    {return _secNode; }
    // ________________________________________________________________________
    void                                   setPrimNode(
        nat p)
    {_primNode = p; }
    // ________________________________________________________________________
    void                                   setSecNode(
        nat p)
    {_secNode = p; }
    // ________________________________________________________________________
    nat                                    getOtherNode(
        nat aNode) const
    {
        assert(hasNode(aNode));
        return aNode == _primNode ? _secNode : _primNode;
    }
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&      s,
        const BranchPlain& c)
    {
        s << c._primNode << "," << c._secNode;
        return s;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    nat _primNode;
    nat _secNode;
};


///////////////////////////////////////////////////////////////////////////////
//                             INLINE DEFINITIONS                            //
///////////////////////////////////////////////////////////////////////////////

namespace std {
template<>
struct hash<BranchPlain>
{
    size_t                                 operator()(
        const BranchPlain& x) const
    {
        return std::hash<size_t>()(x.getSecNode()) ^ std::hash<size_t>()(
            x.getPrimNode());
    }
};


template<>
struct equal_to<BranchPlain>
{
    bool                                   operator()(
        const BranchPlain&a,
        const BranchPlain&b)  const
    {
        return a == b || a.getInverted() == b;
    }
};
}

nat                                        getCommonNode(
    const BranchPlain&oneBranch,
    const BranchPlain&otherBranch);
bool                                       isAdjacent(
    const BranchPlain& oneBranch,
    const BranchPlain& otherBranch);

#endif /* BRANCHBASE_H */

