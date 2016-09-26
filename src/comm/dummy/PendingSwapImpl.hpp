#ifndef _PENDING_SWAP_IMPL_HPP
#define _PENDING_SWAP_IMPL_HPP

#include "AbstractPendingSwap.hpp"

#include <string>
#include <vector>

class ParallelSetup;

///////////////////////////////////////////////////////////////////////////////
//                             PENDING SWAP IMPL                             //
///////////////////////////////////////////////////////////////////////////////
class PendingSwap::Impl : public AbstractPendingSwap
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    Impl(
        SwapElem elem);
    // ________________________________________________________________________
    virtual ~Impl(){}
    // ________________________________________________________________________
    std::vector<char>                    getRemoteData() const;
    // ________________________________________________________________________
    void                                 initialize(
        ParallelSetup&   pl,
        std::vector<char>myChainSer,
        nat              runid);
    // ________________________________________________________________________
    bool                                 allHaveReceived(
        ParallelSetup& pl);
    // ________________________________________________________________________
    bool                                 isFinished();
};


#endif
