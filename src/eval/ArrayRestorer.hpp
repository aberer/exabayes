/**
 * @file ArrayRestorer.hpp
 *
 * @brief allows to restore the lnl array state before evaluating a proposal
 *
 * Note that this doubles the overall memory consumption
 */

#ifndef _LNL_RESTORER_H
#define   _LNL_RESTORER_H

#include "pll.h"
#include "TreeAln.hpp"
#include "PartitionLikelihood.hpp"
#include "ArrayOrientation.hpp"


///////////////////////////////////////////////////////////////////////////////
//                               ARRAY RESTORER                              //
///////////////////////////////////////////////////////////////////////////////
class ArrayRestorer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ArrayRestorer(
        const TreeAln& traln,
        bool           cacheTipTip,
        bool           cacheTipInner);
    // ________________________________________________________________________
    /**
     *  @brief resets all information about switched arrays
     */
    void                                       resetRestorer(
        const TreeAln&   traln,
        ArrayOrientation&curOrient);
    // ________________________________________________________________________
    /**
     *  @brief restore arrays for the given list of partitions
     */
    void                                       restoreSomePartitions(
        TreeAln&                traln,
        const std::vector<bool>&partitions,
        ArrayOrientation&       evalOrientation,
        ArrayReservoir&         res);
    // ________________________________________________________________________
    /**
     *  @brief traverses a subtree and switches arrays, where necessary
     *  @param virtualRoot the root of the subtree
     */
    void                                       traverseAndCache(
        TreeAln&          traln,
        nodeptr           virtualRoot,
        nat               model,
        ArrayOrientation& arrayOrientation,
        ArrayReservoir&   reservoir);
    // ________________________________________________________________________
    void                                       cache(
        TreeAln&               traln,
        nat                    nodeNumber,
        nat                    partitionId,
        const ArrayOrientation&curOrient);
    // ________________________________________________________________________
    void                                       destroyAndForget(
        TreeAln&traln,
        nat     nodeNumber,
        nat     partitionId);
    // ________________________________________________________________________
    void                                       uncache(
        TreeAln&         traln,
        nat              nodeNumber,
        nat              partitionId,
        ArrayOrientation&curOrient,
        ArrayReservoir&  res);
    // ________________________________________________________________________
    void                                       clearMemory(
        ArrayReservoir&res);
    // ________________________________________________________________________
    /**
     *  @brief for the sev-memory saving technique, we also need to remember
     */
    void                                       enableRestoreGapVector()
    {restoresGapVector = true; }
    // ________________________________________________________________________
    void                                       recycleArray(
        TreeAln&        traln,
        nat             nodeNumber,
        nat             partitionId,
        ArrayReservoir& res);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    std::pair<double*, nat>                    removeArray(
        TreeAln&traln,
        nat     num,
        nat     partition);
    // ________________________________________________________________________
    void                                       insertArray(
        TreeAln&traln,
        nat num,
        nat partition,
        std::pair<double*, nat>toInsert);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<PartitionLikelihood> partitionLikelihoods;
    bool                             restoresGapVector;
    ArrayOrientation                 arrayOrientation;
    bool                             cacheTipTip;
    bool                             cacheTipInner;
};


#endif
