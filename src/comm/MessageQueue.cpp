#include "comm/MessageQueue.hpp"

MessageQueue::MessageQueue(const MessageQueue& rhs) 
    : MessageQueue(rhs, std::lock_guard<std::mutex>(rhs._mtx))
  {
  }



MessageQueue::MessageQueue(const MessageQueue &rhs, const std::lock_guard<std::mutex> & lockedMtx)
    : _asyncMessages(rhs._asyncMessages)
{
}




bool MessageQueue::checkMessage(int tag) const 
{
  auto &&lock = std::lock_guard<std::mutex>(_mtx); 
  return _asyncMessages.find(tag) != end( _asyncMessages); 
} 
