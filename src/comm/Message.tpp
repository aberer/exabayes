#include "GlobalVariables.hpp"


template<typename T>
void Message::post(const std::vector<T> &array)
{
  auto &&lock = std::lock_guard<std::mutex>(_mtx); 
  _message.resize(array.size() * sizeof(T));
  std::copy((byte*)array.data(), (byte*) array.data() + array.size() * sizeof(T) , _message.data() );
  _numRead = 0; 
  incrementGeneration(false);
}


template<typename T>
size_t Message::readMessage(T* ptr, size_t offset, size_t length )
{
  auto &&lock = std::lock_guard<std::mutex> (_mtx); 
  assert( (offset + length)  * sizeof(T) <= _message.size() ); 
  std::copy(_message.data() + offset * sizeof(T), _message.data() + (offset + length) * sizeof(T) , (byte*) ptr); 
  ++_numRead; 

  assert(_message.size() % sizeof(T) == 0 ); 
  
  return _message.size() / sizeof(T); 
}


template<typename T>
size_t Message::getMessageSize()  const 
{
  auto &&lock = std::lock_guard<std::mutex> (_mtx);   
  assert(_message.size() % sizeof(T) == 0); 
  return _message.size() / sizeof(T); 
}


template<typename T>
std::vector<T> Message::getMessage()
{
  auto siz = getMessageSize<T>(); 
  auto result = std::vector<T>(siz); 
  auto read = readMessage(result.data(), 0, siz); 
  assert(read == siz); 
  return result; 
} 

