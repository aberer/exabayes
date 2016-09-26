/**
 *  @file Path.hpp represents a  path in a tree
 *
 *  Currently we are assuming, the path always is connected.
 *
 *  This is an ancient class; not entirely happy with it.
 *
 */


#ifndef _PATH_H
#define _PATH_H

#include <vector>

#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "PriorBelief.hpp"

///////////////////////////////////////////////////////////////////////////////
//                                    PATH                                   //
///////////////////////////////////////////////////////////////////////////////
class Path
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    Path()
        : stack(std::vector<BranchPlain>{}){}
    // ________________________________________________________________________
    virtual ~Path(){}

    // ________________________________________________________________________
    /** @brief returns true, if the node with a given id is part of this branch
     * */
    bool                                   nodeIsOnPath(
        int node) const;
    // ________________________________________________________________________
    /** @brief asserts that this path exists in a given tree */
    void                                   debug_assertPathExists(
        TreeAln& traln);
    // ________________________________________________________________________
    /** @brief only add a branch to the path, if it is novel. If the new
     *  branch cancels out an existing branch, the path is shortened again */
    void                                   pushToStackIfNovel(
        BranchPlain   b,
        const TreeAln&traln);
    // ________________________________________________________________________
    // straight-forward container methods
    void                                   append(
        BranchPlain value);
    // ________________________________________________________________________
    void                                   clear();
    // ________________________________________________________________________
    /** @brief number of branches in the path */
    size_t                                 size() const
    {return stack.size(); }
    // ________________________________________________________________________
    /** @brief yields the branch */
    BranchPlain&                           at(
        size_t num){return stack[num]; }
    // ________________________________________________________________________
    BranchPlain                            at(
        size_t num) const {return stack[num]; }
    // ________________________________________________________________________
    /** @brief reverse the path */
    void                                   reverse();
    // ________________________________________________________________________
    /** @brief removes the last element */
    void                                   pop();
    // ________________________________________________________________________
    /** @brief removes the first element */
    void                                   popFront();
    // ________________________________________________________________________
    /** @brief returns the id of the nth node in the path. nodes 0 and n+1 are
     * the outer nodes in this path that do not have a neighbor. */
    nat                                    getNthNodeInPath(
        size_t num) const;
    // ________________________________________________________________________
    /** @brief gets the number of nodes represented by the path (assuming it is
     * connected)  */
    size_t                                 getNumberOfNodes() const
    {return stack.size()  + 1; }
    // ________________________________________________________________________
    void                                   findPath(
        const TreeAln& traln,
        nodeptr        p,
        nodeptr        q);
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&out,
        const Path&  rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // METHODS
    bool                                   findPathHelper(
        const TreeAln&    traln,
        nodeptr           p,
        const BranchPlain&target);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::vector<BranchPlain> stack;
};


#endif
