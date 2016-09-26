#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H

#include "TreeAln.hpp"

class Randomness;
class ParallelSetup;

///////////////////////////////////////////////////////////////////////////////
//                              TREE RANDOMIZER                              //
///////////////////////////////////////////////////////////////////////////////
class TreeRandomizer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    // creates a random tree
    //
    static void                           randomizeTree(
        TreeAln&   traln,
        Randomness&rand);
    // ________________________________________________________________________
    // draws a branch with uniform probability.
    // we have to treat inner and outer branches separatedly.
    //
    static BranchPlain                    drawBranchUniform(
        const TreeAln& traln,
        Randomness&    rand);
    // ________________________________________________________________________
    // samples an inner branch (including orientation), such that each
    // oriented inner branch is equally likely.
    //
    static BranchPlain                    drawInnerBranchUniform(
        const TreeAln& traln,
        Randomness&    rand);
    // ________________________________________________________________________
    static BranchPlain                    drawBranchWithInnerNode(
        const TreeAln& traln,
        Randomness&    rand);
    // ________________________________________________________________________
    static nat                            drawInnerNode(
        const TreeAln& traln,
        Randomness&    rand);

    // ________________________________________________________________________
    // creates a parsimony tree and applies it to the specified tree.
    // the routine does not enforce treatment of branch lengths.
    //
    static void                           createParsimonyTree(
        TreeAln&       traln,
        Randomness&    rand,
        ParallelSetup& pl);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    // currently just exists for the tree randomizer
    //
    // Notice that an inner branch has 3 nodes associated with it,
    // thus the probability of drawing a tip should be accordingly
    // lower.
    //
    static BranchPlain                    drawBranchUniform_helper(
        const TreeAln&traln,
        Randomness&   rand,
        nat           curNumTax);
    // ________________________________________________________________________
    static void                           createStepwiseAdditionParsimonyTree(
        TreeAln&       traln,
        ParallelSetup& pl);
    // ________________________________________________________________________
    static void                           stepwiseAddition(
        pllInstance*   tr,
        partitionList* pr,
        nodeptr        p,
        nodeptr        q,
        ParallelSetup& pl);
};


#endif

