#ifndef _ABSTRACT_PENDING_SWAP_HPP
#define _ABSTRACT_PENDING_SWAP_HPP

#include "SwapElem.hpp"
#include "CommRequest.hpp"
#include "PendingSwap.hpp"

#include <cstdint>

///////////////////////////////////////////////////////////////////////////////
//                           ABSTRACT PENDING SWAP                           //
///////////////////////////////////////////////////////////////////////////////
class AbstractPendingSwap
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    AbstractPendingSwap(
        SwapElem swap)
        : _swap(swap)
    {}
    // ________________________________________________________________________
    virtual ~AbstractPendingSwap()
    {}

    // ________________________________________________________________________
    SwapElem                                     getSwap() const
    {return _swap; }
    // ________________________________________________________________________
    virtual std::vector<char>                    getRemoteData() const = 0;
    // ________________________________________________________________________
    virtual bool                                 isFinished() = 0;
    // ________________________________________________________________________
    virtual bool                                 allHaveReceived(
        ParallelSetup& pl)  = 0;
    // ________________________________________________________________________
    virtual void                                 initialize(
        ParallelSetup&   pl,
        std::vector<char>myChainSer,
        nat              runid)  = 0;
    // ________________________________________________________________________
    static uint64_t                              cantorPair(
        uint64_t a,
        uint64_t b)
    {return (a + b) * (a + b + 1) / 2 + b; }

    ///////////////////////////////////////////////////////////////////////////
    //                            PROTECTED DATA                             //
    ///////////////////////////////////////////////////////////////////////////
protected:
    SwapElem _swap;
};


#endif
