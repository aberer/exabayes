#ifndef _MESSAGE_QUEUE_SINGLE_HPP
#define _MESSAGE_QUEUE_SINGLE_HPP

/**
 *  @brief
 *  this is a true lockless single producer, single consumer message queue.
 *
 *  see sutters article in drdobbs
 */

#include <vector>
#include <limits>
#include <atomic>

using byte =  uint8_t;

///////////////////////////////////////////////////////////////////////////////
//                            MESSAGE QUEUE SINGLE                           //
///////////////////////////////////////////////////////////////////////////////
class MessageQueueSingle
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    struct Node;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    MessageQueueSingle();
    // ________________________________________________________________________
    MessageQueueSingle(
        const MessageQueueSingle& rhs);
    // ________________________________________________________________________
    MessageQueueSingle(
        MessageQueueSingle&& rhs);
    // ________________________________________________________________________
    MessageQueueSingle&                                  operator=(
        MessageQueueSingle rhs);
    // ________________________________________________________________________
    friend void                                          swap(
        MessageQueueSingle& lhs,
        MessageQueueSingle& rhs);
    // ________________________________________________________________________
    ~MessageQueueSingle();
    // ________________________________________________________________________
    template<typename T>
    void                                                 produce(
        std::vector<T>msg);
    // ________________________________________________________________________
    template<typename T>
    std::tuple<bool, std::vector<T> >                    consume(
        int myId);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    Node*              _first; // producer only
    // char pad[CACHE_LINE_SIZE - sizeof(Node*)];
    std::atomic<Node*> _divider; // shared
    // char pad2[CACHE_LINE_SIZE - sizeof(std::atomic<Node*>)];
    std::atomic<Node*> _last;   // shared
};


///////////////////////////////////////////////////////////////////////////////
//                         MESSAGE QUEUE SINGLE NODE                         //
///////////////////////////////////////////////////////////////////////////////
struct MessageQueueSingle::Node
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
    // ________________________________________________________________________
    Node(
        std::vector<byte>msg);
    // ________________________________________________________________________
    Node(
        const Node& rhs) = default;
    // ________________________________________________________________________
    Node&                                                operator=(
        const Node&rhs)  = default;

    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    std::vector<byte>               _message;
    Node*                           _next;
};


#include "MessageQueueSingle.tpp"


#endif
