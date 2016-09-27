#ifndef _SYNC_OSTREAM_HPP
#define _SYNC_OSTREAM_HPP

#include "GlobalVariables.hpp"

#include <sstream>
#include <mutex>

/**
 *  thread safe stream
 */


///////////////////////////////////////////////////////////////////////////////
//                                  SYNC OUT                                 //
///////////////////////////////////////////////////////////////////////////////
class SyncOut
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    SyncOut()
        : _lock(mtx)
        , _ss{}
    {}
    // ________________________________________________________________________
    ~SyncOut(){}
    // ________________________________________________________________________
    template<typename T>
    friend SyncOut&                        operator<<(
        SyncOut& out,
        T        elem)
    {
        out._ss << elem;
        return out;
    }
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&  out,
        const SyncOut& rhs)
    {
        out << rhs._ss.str();
        return out;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::lock_guard<std::mutex>_lock;
    std::stringstream _ss;
};


#endif
