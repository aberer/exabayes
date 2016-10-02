#ifndef FULL_CACHE_POLICY_HPP
#define FULL_CACHE_POLICY_HPP

#include "ArrayPolicy.hpp"
#include "ArrayOrientation.hpp"
#include "ArrayRestorer.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             FULL CACHE POLICY                             //
///////////////////////////////////////////////////////////////////////////////
class FullCachePolicy : public ArrayPolicy
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    FullCachePolicy(
        const TreeAln& traln,
        bool           cacheTipTip,
        bool           cacheTipInner);
    // ________________________________________________________________________
    virtual ~FullCachePolicy(){}
    // ________________________________________________________________________
    virtual void                    imprintPolicy(
        const TreeAln&   traln,
        ArrayOrientation&arrayOrient);
    // ________________________________________________________________________
    virtual void                    freeMemory(
        ArrayReservoir&res);
    // ________________________________________________________________________
    virtual void                    accountForRejectionPolicy(
        TreeAln&                traln,
        const std::vector<bool>&partitions,
        const std::vector<nat>& invalidNodes,
        ArrayOrientation&       arrayOrient,
        ArrayReservoir&         res);
    // ________________________________________________________________________
    virtual void                    prepareForEvaluation(
        TreeAln&          traln,
        BranchPlain       virtualRoot,
        nat               model,
        ArrayOrientation& arrayOrientation,
        ArrayReservoir&   res);
    // ________________________________________________________________________
    virtual auto                    clone() const
        ->ArrayPolicy::UPtr;
    // ________________________________________________________________________
    void                            enableRestoreGapVector()
    {restorer.enableRestoreGapVector();}

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    ArrayRestorer restorer;
};


#endif
