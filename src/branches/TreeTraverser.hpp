#ifndef _TREE_TRAVERSER_HPP
#define _TREE_TRAVERSER_HPP

#include "LikelihoodEvaluator.hpp"
#include "TreeAln.hpp"
#include "InsertionResult.hpp"
#include "ParallelSetup.hpp"
#include "BranchLengthOptimizer.hpp"

#include <functional>
#include <unordered_map>

///////////////////////////////////////////////////////////////////////////////
//                                 EVAL PLAIN                                //
///////////////////////////////////////////////////////////////////////////////
class EvalPlain
{
public:
    // ________________________________________________________________________
    InsertionResult                                 getResult(
        LikelihoodEvaluator&                  eval,
        TreeAln&                              traln,
        const BranchPlain&                    insertBranch,
        const BranchPlain&                    prunedSubtree,
        const BranchPlain&                    draggedBranch,
        const std::vector<AbstractParameter*>&params) const;
};


///////////////////////////////////////////////////////////////////////////////
//                              EVAL OPT DRAGGED                             //
///////////////////////////////////////////////////////////////////////////////
class EvalOptDragged
{
public:
    // ________________________________________________________________________
    InsertionResult                                 getResult(
        LikelihoodEvaluator&                  eval,
        TreeAln&                              traln,
        const BranchPlain&                    insertBranch,
        const BranchPlain&                    prunedSubtree,
        const BranchPlain&                    draggedBranch,
        const std::vector<AbstractParameter*>&params) const;
};


///////////////////////////////////////////////////////////////////////////////
//                               EVAL OPT THREE                              //
///////////////////////////////////////////////////////////////////////////////
class EvalOptThree
{
public:
    // ________________________________________________________________________
    InsertionResult                                 getResult(
        LikelihoodEvaluator&                  eval,
        TreeAln&                              traln,
        const BranchPlain&                    insertBranch,
        const BranchPlain&                    prunedSubtree,
        const BranchPlain&                    draggedBranch,
        const std::vector<AbstractParameter*>&params) const;
};


///////////////////////////////////////////////////////////////////////////////
//                               TREE TRAVERSER                              //
///////////////////////////////////////////////////////////////////////////////
template<class InsertEval>
class TreeTraverser
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    TreeTraverser(
        int                            depth,
        TreeAln&                       traln,
        LikelihoodEvaluator&           eval,
        std::vector<AbstractParameter*>params,
        const BranchPlain&             rootOfTraversed,
        const BranchPlain&             prunedSubtree);
    // ________________________________________________________________________
    std::vector<InsertionResult>                    getResult() const
    {return _result; }
    // ________________________________________________________________________
    void                                            traverse();

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                                            testInsert(
        const BranchPlain&insertBranch,
        BranchLengths     floatingBranch,
        bool              isFirst,
        int               curDepth);
    ///////////////////////////////////////////////////////////////////////////
    //                             PRIVATE DATA                              //
    ///////////////////////////////////////////////////////////////////////////
private:
    bool                                        _doFirst;
    int                                         _depth;
    std::reference_wrapper<TreeAln>             _traln;
    std::reference_wrapper<LikelihoodEvaluator> _eval;
    std::vector<AbstractParameter*>             _params;
    BranchPlain                                 _rootOfTraversed;
    BranchPlain                                 _prunedSubtree;
    std::vector<InsertionResult>                _result;

    static const int                            _maxIter;

    InsertEval                                  _insertEval;
};


#include  "TreeTraverser.tpp"

#endif
