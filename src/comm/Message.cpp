#include "comm/Message.hpp"
#include "GlobalVariables.hpp"



Message::Message()
  : _generation{0}
  , _numRead{0}
  {
  }


Message::Message(const Message& rhs, const std::lock_guard<std::mutex> &lock)
  : _generation{rhs._generation}
  , _numRead{rhs._numRead}
{
  
}


Message::Message(const Message& rhs) 
  : Message(rhs, std::lock_guard<std::mutex>(rhs._mtx))
{
}


bool Message::clearMessage(int numToBeRead)
{
  auto&& lock = std::lock_guard<std::mutex>(_mtx); 

  if(_numRead == numToBeRead)
    {
      _numRead = 0; 
      return true; 
    }
      
  return false; 
}


void Message::incrementGeneration(bool lock)
{
  if(lock )
    _mtx.lock(); 

  if(_generation == std::numeric_limits<uint8_t>::max())
    _generation = 0; 
  else 
    ++_generation; 

  if(lock)
    _mtx.unlock();
}


uint8_t Message::getGeneration() const 
{
  auto &&lock = std::lock_guard<std::mutex>(_mtx); 
  return _generation; 
}



