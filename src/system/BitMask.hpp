#ifndef _BITMASK_HPP
#define _BITMASK_HPP


template<typename T> 
class BitMask
{
public: 
  static int constexpr length = sizeof(T) * 8 ; // dont know where size_of_byte define is 

public: 
  BitMask(); 
  T operator[](int num ) {return _impl[num]; } 
private:
  T _impl[length]; 
}; 


#endif
