#ifndef CONSUMABLE_MESSAGE_HPP
#define CONSUMABLE_MESSAGE_HPP

#include <vector>
#include <mutex>

#include <limits>

typedef uint8_t byte; 

class ConsumableMessage
{
public: 
  template<typename T >
  ConsumableMessage(std::vector<T> data, int numReads); 
  
  ConsumableMessage (const ConsumableMessage& rhs)
    : ConsumableMessage(rhs, std::lock_guard<std::mutex>(rhs._mtx))
  {
  } 

  
  // template<typename T>
  // void initialize(std::vector<T> data, int numReads); 

  template<typename T>
  std::tuple<std::vector<T>,int> consume() ; 

private : 			// METHODS 
  ConsumableMessage( const ConsumableMessage& rhs,  const std::lock_guard<std::mutex> &lockedMtx)
    : _message(rhs._message)
    , _readsLeft(rhs._readsLeft)
  {
  } 

private: 			// ATTRIBUTES
  std::vector<byte> _message; 
  int _readsLeft; 
  mutable std::mutex _mtx; 
}; 

#include "comm/ConsumableMessage.tpp"



#endif
