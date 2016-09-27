#ifndef _MOVE_ENUMERATOR_HPP
#define _MOVE_ENUMERATOR_HPP

#include "common.h"

#include "TreeAln.hpp"

#include <vector>

class SprMove;

///////////////////////////////////////////////////////////////////////////////
//                              MOVE ENUMERATION                             //
///////////////////////////////////////////////////////////////////////////////
class MoveEnumeration
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    static auto                    getAllUniqueMoves(
        const TreeAln& traln,
        nat            dist)
        ->std::vector<SprMove>;
    // ________________________________________________________________________
    static auto                    getMovesForBranch(
        const TreeAln& traln,
        BranchPlain    pruneBranch,
        nat            dist)
        ->std::vector<SprMove>;
    // ________________________________________________________________________
    static auto                    getMovesFromDepthFirstSearch(
        const TreeAln&     traln,
        const BranchPlain& pruneBranch,
        int                dist)
        ->std::vector<SprMove>;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    static auto                    getMovesFromDepthFirstSearch_helper(
        const TreeAln&     traln,
        const BranchPlain& pruneBranch,
        const BranchPlain& insert,
        int                depth,
        bool               isFirst)
        ->std::vector<SprMove>;
};


#endif
