#ifndef _NEW_COMMUNICATOR_HPP
#define _NEW_COMMUNICATOR_HPP

#include "RemoteComm.hpp"
#include "LocalComm.hpp"
#include "threadDefs.hpp"

#include <vector>
#include <cassert>

#include <unordered_map>

class ThreadResource;

///////////////////////////////////////////////////////////////////////////////
//                                COMMUNICATOR                               //
///////////////////////////////////////////////////////////////////////////////
class Communicator
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
    using SELF =  Communicator;
public:
    // ________________________________________________________________________
    /**
     *  @param tid2rank absolute ranks in the communicator ; contains only
     * ranks of local threads
     */
    Communicator(
        std::unordered_map<tid_t, int>tid2rank);
    // ________________________________________________________________________
    Communicator(
        const Communicator& rhs) = delete;
    // ________________________________________________________________________
    Communicator(
        Communicator&& rhs) = default;
    // ________________________________________________________________________
    ~Communicator(){}
    // ________________________________________________________________________
    Communicator&                          operator=(
        Communicator rhs);
    // ________________________________________________________________________
    friend void                            swap(
        Communicator& lhs,
        Communicator& rhs);
    // ________________________________________________________________________
    void                                   createSendRequest(
        std::vector<char>array,
        int              dest,
        int              tag,
        CommRequest&     req);
    // ________________________________________________________________________
    void                                   createRecvRequest(
        int          src,
        int          tag,
        nat          length,
        CommRequest& req);

#include "CommCore.hpp"

    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&       out,
        const Communicator& rhs);
    // ________________________________________________________________________
    int                                    mapToLocalRank(
        int rank) const;
    // ________________________________________________________________________
    int                                    mapToRemoteRank(
        int rank) const;
    // ________________________________________________________________________
    int                                    getProcsPerNode();
    // ________________________________________________________________________
    LocalComm&                             getLocalComm();
    // ________________________________________________________________________
    RemoteComm&                            getRemoteComm();
    // ________________________________________________________________________
    void                                   initWithMaxChains(
        size_t numChains,
        size_t numThreadsChecking);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    RemoteComm                             _remoteComm;
    LocalComm                              _localComm;
};


#include "Communicator.tpp"

#endif

