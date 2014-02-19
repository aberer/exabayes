#ifndef _MESSAGE_QUEUE 
#define _MESSAGE_QUEUE 

#include "comm/ConsumableMessage.hpp"
#include <unordered_map>

class MessageQueue
{
public: 
  MessageQueue(){}

  MessageQueue(const MessageQueue& rhs); 

  template<typename T>
  void post(const std::vector<T> &message, int numRead, int tag);  

  template<typename T>
  std::tuple<bool,std::vector<T>> consumeOne(int tag); 

  bool checkMessage(int tag) const ; 


private: 			// METHODS 
  MessageQueue(const MessageQueue &rhs, const std::lock_guard<std::mutex> & lockedMtx); 

private: 			// ATTRIBUTES
  std::unordered_map<int,ConsumableMessage> _asyncMessages; 
  mutable std::mutex _mtx; 
}; 


#include "comm/MessageQueue.tpp"

#endif
