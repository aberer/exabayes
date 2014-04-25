#ifndef _LOCK_FREE_QUEUE_HPP
#define _LOCK_FREE_QUEUE_HPP


/** 
    Notice: this data structure is not truely lock-free...

    because of false-sharing, this structure may even be faster, if
    multiple messages are produced (one per consumer)

    it is a modified version of 
    http://www.drdobbs.com/parallel/writing-a-generalized-concurrent-queue/211601363?pgno=4
    that mostly avoids some copies 
 */ 




#include <vector>
#include <limits>
#include <atomic>


typedef uint8_t byte;  

class  MessageQueue 
{
  struct Node; 

public: 
  explicit MessageQueue(int numThreads); 
  MessageQueue(const MessageQueue& rhs)  ; 
  MessageQueue( MessageQueue&& rhs)   ; 
  MessageQueue& operator=(MessageQueue rhs); 
  friend void swap(MessageQueue& lhs, MessageQueue& rhs); 
  ~MessageQueue(); 
  
  template<typename T>
  void produce(std::vector<T> msg, const std::vector<int> &forWhom); 
  template<typename T>
  std::tuple<bool,std::vector<T> >  consume(int myId); 

  
private: 
  Node* _first;			// producer only 
  // char pad[CACHE_LINE_SIZE - sizeof(Node*)]; 
  std::atomic<Node*> _divider; 	// shared 
  // char pad2[CACHE_LINE_SIZE - sizeof(std::atomic<Node*>)]; 
  std::atomic<Node*> _last; 	// shared 
  // char pad2[CACHE_LINE_SIZE - sizeof(std::atomic<Node*>)]; 
  std::atomic<bool> _consumerLock;
  // char pad2[CACHE_LINE_SIZE - sizeof(std::atomic<bool>)]; 
  int _numThreads; 
}; 

struct MessageQueue::Node {
  // METHODS
  Node(std::vector<byte> msg, std::vector<int> whoReads ); 
  std::tuple<int, std::vector<byte> > getConsumed(int threadId); 

  // ATTRIBUTES
  std::vector<byte> _message; 
  std::atomic<int> _numConsume; 
  Node *_next; 
  std::vector<int>  _whoReads; 
}; 


#include "comm/threads/MessageQueue.tpp"


#endif
