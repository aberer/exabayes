#ifndef NO_CACHE_POLICY_HPP
#define NO_CACHE_POLICY_HPP

#include "ArrayPolicy.hpp"
#include "BranchPlain.hpp"

///////////////////////////////////////////////////////////////////////////////
//                              NO CACHE POLICY                              //
///////////////////////////////////////////////////////////////////////////////
class NoCachePolicy : public ArrayPolicy
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    NoCachePolicy(
        const TreeAln&traln);
    // ________________________________________________________________________
    virtual ~NoCachePolicy(){}
    // ________________________________________________________________________
    virtual void
                        imprintPolicy(
        const TreeAln&   traln,
        ArrayOrientation&arrayOrient){}
    // ________________________________________________________________________
    virtual void
                        freeMemory(
        ArrayReservoir& res){}
    // ________________________________________________________________________
    virtual void
                        prepareForEvaluation(
        TreeAln&          traln,
        BranchPlain       virtualRoot,
        nat               models,
        ArrayOrientation& arrayOrientation,
        ArrayReservoir&   res)
    {}
    // ________________________________________________________________________
    virtual auto
                        accountForRejectionPolicy(
        TreeAln&                traln,
        const std::vector<bool>&partitions,
        const std::vector<nat>& invalidNodes,
        ArrayOrientation&       arrayOrient,
        ArrayReservoir&         res)
        ->void;
    // ________________________________________________________________________
    virtual std::unique_ptr<ArrayPolicy>
                        clone() const;
};


#endif
