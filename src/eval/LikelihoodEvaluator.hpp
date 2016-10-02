#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "ArrayPolicy.hpp"
#include "TreeAln.hpp"
#include "ArrayOrientation.hpp"

#include "extensions.hpp"

#include <vector>

class ParallelSetup;

///////////////////////////////////////////////////////////////////////////////
//                            LIKELIHOOD EVALUATOR                           //
///////////////////////////////////////////////////////////////////////////////
class LikelihoodEvaluator
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    LikelihoodEvaluator(
        const TreeAln&                 traln,
        ArrayPolicy*                   plcy,
        std::shared_ptr<ArrayReservoir>arrayReservoir,
        ParallelSetup*                 pl);
    // ________________________________________________________________________
    ~LikelihoodEvaluator(){}
    // ________________________________________________________________________
    LikelihoodEvaluator(
        LikelihoodEvaluator&&rhs) = default;
    // ________________________________________________________________________
    LikelihoodEvaluator(
        const LikelihoodEvaluator&rhs);
    // ________________________________________________________________________
    LikelihoodEvaluator&                   operator=(
        LikelihoodEvaluator rhs);
    // ________________________________________________________________________
    friend void                            swap(
        LikelihoodEvaluator&lhs,
        LikelihoodEvaluator&rhs);
    // ________________________________________________________________________
    /**
     * @brief: evaluate a list of partitions. This is always a full traversal
     */
    void                                   evaluatePartitionsWithRoot(
        TreeAln&                traln,
        const BranchPlain&      root,
        const std::vector<nat>& partitions,
        bool                    fullTraversal);
    // ________________________________________________________________________
    /**
     *  @brief a horrible hack: since I do not get what arcane
     *  invocations are necessary to make branch length optimization
     *  work without an evaluate, we execute an evaluate here, but with
     *  partition.width set to 0. Computationally this is okay, but
     *  otherwise it is horrible. If you have a few days to spent,
     *  please correct this.
     */
    void                                   evaluatePartitionsDry(
        TreeAln&                traln,
        const BranchPlain&      root,
        const std::vector<nat>& partitions);
    // ________________________________________________________________________
    void                                   invalidateWithBranches(
        const TreeAln&                        traln,
        const BranchPlain&                    root,
        const std::unordered_set<BranchPlain>&invalidBranches);
    // ________________________________________________________________________
    /**
     *  @brief evaluate a subtree. This used to be the
     *  newview-command. The this-node is the important part, the
     *  that-node of the branch only specifies the orientation.
     */
    void                                   evalSubtree(
        TreeAln&          traln,
        nat               partition,
        const BranchPlain&evalBranch,
        bool              fullTraversal);
    // ________________________________________________________________________
    /**
     *  @brief evaluation at a given branch
     */
    void                                   evaluate(
        TreeAln&          traln,
        const BranchPlain&root,
        bool              fullTraversal);
    // ________________________________________________________________________
    /**
     *  @brief mark a node as dirty
     */
    void                                   invalidateArray(
        const TreeAln&traln,
        nat           nodeId);
    // ________________________________________________________________________
    void                                   invalidateArray(
        const TreeAln&traln,
        nat           partitionId,
        nat           nodeId);
    // ________________________________________________________________________
    void                                   debugPrintToCompute(
        const TreeAln&    traln,
        const BranchPlain&root);
    // ________________________________________________________________________
    void                                   debugPrintToComputeHelper(
        const TreeAln&    traln,
        const BranchPlain&root);
    // ________________________________________________________________________
    /**
     *  @brief marks the entire tree dirty
     */
    void                                   markPartitionDirty(
        const TreeAln&traln,
        nat           partition);
    // ________________________________________________________________________
    /**
     *  @brief indicates whether a likelihood array for a given node and
     *  partition is invalid
     */
    bool                                   isDirty(
        nat partition,
        nat nodeId)  const;
    // ________________________________________________________________________
    /**
     *  @brief make the current state in the tree resettable (only needed for
     * chain.cpp)
     */
    void                                   imprint(
        const TreeAln&traln);
    // ________________________________________________________________________
    /**
     *  @brief the main point of this is to free remaining memory
     */
    void                                   freeMemory();
    // ________________________________________________________________________
    /**
     *  @brief use backup likelihood
     */
    void                                   accountForRejection(
        TreeAln&                traln,
        const std::vector<bool>&partitions,
        const std::vector<nat>& invalidNodes);
    // ________________________________________________________________________
    ArrayOrientation                       getOrientation() const
    {return _arrayOrientation;}
    // ________________________________________________________________________
    void                                   expensiveVerify(
        TreeAln&    traln,
        BranchPlain root,
        log_double  toVerify);
    // ________________________________________________________________________
    void                                   setDebugTraln(
        std::shared_ptr<TreeAln>_debugTraln);
    // ________________________________________________________________________
    ArrayReservoir&                        getArrayReservoir()
    {return *(_arrayReservoir);}
    // ________________________________________________________________________
    ParallelSetup&                         getParallelSetup()
    {return *_plPtr;}
    // ________________________________________________________________________
    void                                   evaluateSubtrees(
        TreeAln&                traln,
        const BranchPlain&      root,
        const std::vector<nat>& partitions,
        bool                    fullTraversal);
    // ________________________________________________________________________
    void                                   evaluateSubtrees(
        TreeAln&          traln,
        const BranchPlain&root,
        bool              fullTraversal);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    /**
     *  @brief traverses the tree and reorients nodes to enforce
     *  re-computation
     *
     *  @param the new virtual root
     *
     *  @notice an important assumption is that the new virtual root is
     *  adjacent to invalidated likelihood arrays. Everything else would
     *  be extremely inefficient.
     */
    bool                                   applyDirtynessToSubtree(
        TreeAln&          traln,
        nat               partition,
        const BranchPlain&branch);
    // ________________________________________________________________________
    void                                   disorientDebug(
        TreeAln&           traln,
        const BranchPlain& root);
    // ________________________________________________________________________
    void                                   disorientDebugHelper(
        TreeAln&           traln,
        const BranchPlain& root);
    // ________________________________________________________________________
    bool                                   subtreeValid(
        const TreeAln&                        traln,
        const BranchPlain&                    root,
        const std::unordered_set<BranchPlain>&invalidBranches);
    // ________________________________________________________________________
    void                                   evaluateLikelihood(
        TreeAln&    traln,
        BranchPlain branch,
        bool        fullTraversal,
        bool        isDebugEval);
    // ________________________________________________________________________
    void                                   updatePartials(
        TreeAln& traln,
        nodeptr  p);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::shared_ptr<TreeAln>_debugTralnPtr;
    bool _verifyLnl;
    ArrayPolicy::UPtr _arrayPolicy;
    ArrayOrientation _arrayOrientation;
    std::shared_ptr<ArrayReservoir>_arrayReservoir;
    ParallelSetup*                  _plPtr;
};


#endif

