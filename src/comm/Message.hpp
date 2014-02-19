#ifndef _MY_MESSAGE_HPP
#define _MY_MESSAGE_HPP

#include "common.h"
#include <iostream>
#include <cassert>
#include <mutex>

#include <cstdint>
#include <algorithm>

typedef  uint8_t byte; 



class Message
{
public: 
  Message(); 
  Message(const Message& rhs) ; 
  ~Message(){}

  template<typename T>
  void post(const std::vector<T> &array) ; 
  
  template<typename T>
  size_t getMessageSize() const ; 

  template<typename T>
  size_t readMessage(T* ptr, size_t offset, size_t length )  ; 
  
  template<typename T>
  std::vector<T> getMessage(); 

  bool clearMessage(int numToBeRead); 
  uint8_t getGeneration() const ; 
  void incrementGeneration(bool lock); 

private: 			// METHODS 
  Message(const Message& rhs, const std::lock_guard<std::mutex> &lock);



private: 			// ATTRIBUTES
  std::vector<byte> _message; 
  uint8_t _generation; 
  int _numRead; 
  mutable std::mutex _mtx; 
}; 


#include "comm/Message.tpp"

#endif
