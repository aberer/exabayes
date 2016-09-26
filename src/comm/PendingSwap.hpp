#ifndef _UNFINISHED_SWAP
#define _UNFINISHED_SWAP

#include <cassert>
#include <vector>
#include <memory>

#include "common.h"

class ParallelSetup;
class SwapElem;
class AbstractPendingSwap;

///////////////////////////////////////////////////////////////////////////////
//                                PENDING SWAP                               //
///////////////////////////////////////////////////////////////////////////////
class PendingSwap
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    class Impl;         // since we have an abstract variant, this does not
                        // really fit anymore
    // ________________________________________________________________________
    PendingSwap(
        SwapElem swap,
        bool     isLocal);
    // ________________________________________________________________________
    PendingSwap(
        PendingSwap&& rhs);
    // ________________________________________________________________________
    PendingSwap&                         operator=(
        const PendingSwap& elem)   = delete;
    // ________________________________________________________________________
    PendingSwap(
        const PendingSwap& rhs) = delete;
    // ________________________________________________________________________
    ~PendingSwap();
    // ________________________________________________________________________
    friend void                          swap(
        PendingSwap& lhs,
        PendingSwap& rhs);
    // ________________________________________________________________________
    bool                                 isFinished();
    // ________________________________________________________________________
    void                                 initialize(
        ParallelSetup&   pl,
        std::vector<char>myChainSer,
        nat              runid);
    // ________________________________________________________________________
    SwapElem                             getSwap() const;
    // ________________________________________________________________________
    bool                                 allHaveReceived(
        ParallelSetup& pl);
    // ________________________________________________________________________
    std::vector<char>                    getRemoteData()  const;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    int                                  createTag(
        nat numChains) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::unique_ptr<AbstractPendingSwap> _impl;
};


#endif

