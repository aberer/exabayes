#ifndef _TEE_STREAM_H
#define _TEE_STREAM_H

#include "teebuf.hpp"

#include <iostream>
#include <unordered_set>
#include <vector>

template<typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> elems)
{
  bool isFirst = true; 
  for(auto &elem : elems)
    {
      if(isFirst) 
	isFirst = false; 
      else 
	out << "," ; 
      out << elem ; 
    }
  return out; 
}


template<typename T>
std::ostream& operator<<(std::ostream& out, std::unordered_set<T> elems)
{
  bool isFirst = true; 
  for(auto &elem : elems)
    {
      if(isFirst)
	isFirst = false; 
      else 
	out << "," ; 
      out << elem; 
    }

  return out; 
}


class teestream : public std::ostream
{
public:   
  teestream(std::ostream &o1, std::ostream &o2)
    : std::ostream(&tbuf)
    ,tbuf(o1.rdbuf(), o2.rdbuf())
  {}
  
  void disable(){tbuf.disable();}
  
private :
  teebuf tbuf; 
}; 

#endif
