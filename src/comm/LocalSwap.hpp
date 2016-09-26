#ifndef _LOCAL_SWAP_HPP
#define _LOCAL_SWAP_HPP

#include "AbstractPendingSwap.hpp"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//                                 LOCAL SWAP                                //
///////////////////////////////////////////////////////////////////////////////
class LocalSwap :  public AbstractPendingSwap
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    LocalSwap(
        SwapElem elem)
        : AbstractPendingSwap(elem)
        , _plPtr{nullptr}
        , _haveReceived(false)
        , _dataReceived{}
        , _runid{0}
    {}
    // ________________________________________________________________________
    LocalSwap&                                   operator=(
        const LocalSwap&rhs);
    // ________________________________________________________________________
    // implements AbstractPendingSwap
    virtual std::vector<char>                    getRemoteData()  const;
    // ________________________________________________________________________
    // implements AbstractPendingSwap
    virtual bool                                 isFinished();
    // ________________________________________________________________________
    // implements AbstractPendingSwap
    virtual bool                                 allHaveReceived(
        ParallelSetup& pl);
    // ________________________________________________________________________
    // implements AbstractPendingSwap
    virtual void                                 initialize(
        ParallelSetup&   pl,
        std::vector<char>myChainSer,
        nat              runid);

    ///////////////////////////////////////////////////////////////////////////
    //                             PRIVATE DATA                              //
    ///////////////////////////////////////////////////////////////////////////
private:
    // slight hack, because i was too lazy to change the signatures...
    ParallelSetup*    _plPtr;

    bool              _haveReceived;
    std::vector<char> _dataReceived;
    nat               _runid;
};


#endif
