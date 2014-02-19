#include <cassert>  
#include <limits>

typedef uint8_t byte ; 

template<typename T>
ConsumableMessage::ConsumableMessage(std::vector<T> data, int numReads)
{
  auto &&lock = std::lock_guard<std::mutex>{_mtx}; 
  _message.resize(sizeof(T) * data.size()); 
  std::copy(reinterpret_cast<byte*>(data.data() ) ,reinterpret_cast<byte*>(data.data()) + sizeof(T) * data.size() , begin(_message)); 
  _readsLeft = numReads; 
}


template<typename T>
std::tuple<std::vector<T>,int> ConsumableMessage::consume() 
{
  assert(_message.size() % sizeof(T) == 0); 
  auto result = std::vector<T>{}; 
  result.resize(_message.size() / sizeof(T));
  std::copy(begin(_message), end(_message), reinterpret_cast<byte*>(result.data())); 

  // lock not necesary before ... 
  auto &&lock = std::lock_guard<std::mutex>{_mtx}; 
  assert(_readsLeft > 0); 
  --_readsLeft; 

  return make_tuple(result,_readsLeft);
}  

