#ifndef COMMCORENEW_H
#define COMMCORENEW_H

#include <cassert>

// this class *could* be used for implementing CRTP with communicators.
// However, the overhead is enormous. Instead use the more clumsy CommCore.hpp


///////////////////////////////////////////////////////////////////////////////
//                               COMM CORE NEW                               //
///////////////////////////////////////////////////////////////////////////////
template<class DERIVED>
class CommCoreNew
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    gather(
        std::vector<T>myData,
        nat           root = 0)
    {return static_cast<DERIVED*>(this)->gather_impl(myData, root); }
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    gather_impl(
        std::vector<T>myData,
        nat           root = 0)
    {assert(0); return {}; }
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    gatherVariableLength(
        std::vector<T>myData,
        int           root = 0)
    {
        return static_cast<DERIVED*>(this)
                   ->gatherVariableLength_impl(myData, root);
    }
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    gatherVariableKnownLength(
        std::vector<T>   myData,
        std::vector<int>&countsPerProc,
        std::vector<int>&displPerProc,
        int              root = 0)
    {
        return static_cast<DERIVED*>(this)->gatherVariableKnownLength_impl(
            myData, countsPerProc, displPerProc, root);
    }
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    scatterVariableKnownLength(
        std::vector<T>   allData,
        std::vector<int>&countsPerProc,
        std::vector<int>&displPerProc,
        int              root)
    {
        return static_cast<DERIVED*>(this)->scatterVariableKnownLength_impl(
            allData, countsPerProc, displPerProc, root);
    }
    // ________________________________________________________________________
    template<typename T>
    // done
    std::vector<T>                    broadcast(
        std::vector<T>array,
        int           root = 0)
    {
        return static_cast<DERIVED*>(this)->broadcast_impl(array, root);
    }
    // ________________________________________________________________________
    template<typename T>
    // ok
    T                                 broadcastVar(
        T   var,
        int root = 0)
    {
        return static_cast<DERIVED*>(this)->broadcastVar_impl(var, root);
    }
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    allReduce(
        std::vector<T>myValues)       //
    {
        return static_cast<DERIVED*>(this)->allReduce_impl(myValues);
    }
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                    reduce(
        std::vector<T>data,
        int           root)
    {
        return static_cast<DERIVED*>(this)->reduce_impl(data, root);
    }
    // ________________________________________________________________________
    template<typename T>
    T                                 receive(
        int source,
        int tag)
    {return static_cast<DERIVED*>(this)->receive_impl(source, tag); }
    // ________________________________________________________________________
    template<typename T>
    void                              send(
        T   elem,
        int dest,
        int tag)
    {static_cast<DERIVED*>(this)->send_impl(elem, dest, tag); }
    // ________________________________________________________________________
    int                               getRank() const
    {return static_cast<DERIVED*>(this)->getRank_impl(); }
    // ________________________________________________________________________
    size_t                            size() const
    {return static_cast<DERIVED*>(this)->size_impl(); }
    // ________________________________________________________________________
    bool                              isValid() const
    {return static_cast<DERIVED*>(this)->isValid_impl(); }
    // ________________________________________________________________________
    bool                              haveThreadSupport() const
    {return static_cast<DERIVED*>(this)->haveThreadSupport_impl(); }
    // ________________________________________________________________________
    DERIVED                           split(
        const std::vector<int>&color,
        const std::vector<int>&rank) const
    {return static_cast<DERIVED*>(this)->split_impl(color, rank); }
    // ________________________________________________________________________
    void                              waitAtBarrier()
    {static_cast<DERIVED*>(this)->waitAtBarrier_impl(); }
    // ________________________________________________________________________
    static void                       finalize()
    {DERIVED::finalize_impl(); }
    // ________________________________________________________________________
    static void                       initComm(
        int   argc,
        char**argv)
    {DERIVED::initComm_impl(argc, argv); }
    // ________________________________________________________________________
    static void                       abort(
        int  code,
        bool waitForAll)
    {DERIVED::abort_impl(code, waitForAll); }
};


#endif /* COMMCORENEW_H */
