#ifndef _ARRAY_ORIENTATION_HPP
#define _ARRAY_ORIENTATION_HPP

#include <vector>
#include <iostream>

#include "TreeAln.hpp"

#define INVALID 0

///////////////////////////////////////////////////////////////////////////////
//                             ARRAY ORIENTATION                             //
///////////////////////////////////////////////////////////////////////////////
class ArrayOrientation
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    // decides how far down the tree we search to determine the
    // correctness of the array that belongs to a moved subtree
    // (this is a hack)
    const static nat maxSprSearch;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ArrayOrientation(
        const TreeAln&traln);
    // ________________________________________________________________________
    bool                                   isCorrect(
        nat part,
        nat id,
        nat value) const
    {return orientation.at(part).at(id) ==  value; }
    // ________________________________________________________________________
    bool                                   searchNodeInSubtree(
        const TreeAln&traln,
        nodeptr       p,
        nat           nodeId,
        nat           depth) const;
    // ________________________________________________________________________
    bool                                   isCorrectNew(
        const TreeAln& traln,
        nat            part,
        nodeptr        p) const;
    // ________________________________________________________________________
    nat                                    getOrientation(
        nat part,
        nat id) const;
    // ________________________________________________________________________
    void                                   setOrientation(
        nat part,
        nat id,
        nat value){orientation[part][id] = value; }
    // ________________________________________________________________________
    void                                   setOrientationForAllPartitions(
        nat id,
        nat value);
    // ________________________________________________________________________
    void                                   setPartitionInvalid(
        nat part);
    // ________________________________________________________________________
    void                                   setInvalid(
        nat part,
        nat id);
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&          out,
        const ArrayOrientation&rhs);
    // ________________________________________________________________________
    static nat                             nodeNum2ArrayNum(
        nat nodeNum,
        nat numTax);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<std::vector<nat> > orientation;
};


#endif
