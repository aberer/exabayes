#ifndef _TEE_STREAM_H
#define _TEE_STREAM_H

#include "teebuf.hpp"

#include <iostream>

using namespace std; 

class teestream : public std::ostream
{
public:   
  teestream(ostream &o1, ostream &o2)
    : ostream(&tbuf)
    ,tbuf(o1.rdbuf(), o2.rdbuf())
  {}
  
  void disable(){tbuf.disable();}
  
private :
  teebuf tbuf; 
}; 

#endif
