#ifndef _TEE_STREAM_H
#define _TEE_STREAM_H

#include "teebuf.hpp"

#include <iostream>


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
