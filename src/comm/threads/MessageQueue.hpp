#ifndef _LOCK_FREE_QUEUE_HPP
#define _LOCK_FREE_QUEUE_HPP


/**
 *  Notice: this data structure is not truely lock-free...
 *
 *  because of false-sharing, this structure may even be faster, if
 *  multiple messages are produced (one per consumer)
 *
 *  it is a modified version of
 *
 * http://www.drdobbs.com/parallel/writing-a-generalized-concurrent-queue/211601363?pgno=4
 *  that mostly avoids some copies
 */


#include <vector>
#include <limits>
#include <atomic>


typedef uint8_t byte;

///////////////////////////////////////////////////////////////////////////////
//                               MESSAGE QUEUE                               //
///////////////////////////////////////////////////////////////////////////////
class MessageQueue
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    class Node;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    explicit MessageQueue(
        int numThreads);
    // ________________________________________________________________________
    MessageQueue(
        const MessageQueue& rhs);
    // ________________________________________________________________________
    MessageQueue(
        MessageQueue&& rhs);
    // ________________________________________________________________________
    MessageQueue&                                          operator=(
        MessageQueue rhs);
    // ________________________________________________________________________
    friend void                                            swap(
        MessageQueue& lhs,
        MessageQueue& rhs);
    // ________________________________________________________________________
    ~MessageQueue();
    // ________________________________________________________________________
    template<typename T>
    void                                                   produce(
        std::vector<T>         msg,
        const std::vector<int>&forWhom);
    // ________________________________________________________________________
    template<typename T>
    std::tuple<bool, std::vector<T> >                      consume(
        int myId);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    Node*              _first; // producer only
    std::atomic<Node*> _divider; // shared
    std::atomic<Node*> _last;   // shared
    std::atomic<bool>  _consumerLock;
    int                _numThreads;
};


///////////////////////////////////////////////////////////////////////////////
//                             MESSAGE QUEUE NODE                            //
///////////////////////////////////////////////////////////////////////////////
class MessageQueue::Node
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    Node(
        std::vector<byte>msg,
        std::vector<int> whoReads);
    // ________________________________________________________________________
    Node(
        const Node&node) = default;
    // ________________________________________________________________________
    Node&                                                  operator=(
        const Node&rhs)  = default;
    // ________________________________________________________________________
    std::tuple<int, std::vector<byte> >                    getConsumed(
        int threadId);


    // TODO: omg
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    std::vector<byte> _message;
    std::atomic<int>  _numConsume;
    Node*             _next;
    std::vector<int>  _whoReads;
};


#include "MessageQueue.tpp"


#endif
