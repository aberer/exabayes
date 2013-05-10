#include "teebuf.hpp"

#include <iostream>

using namespace std; 

class teestream : public std::ostream
{
public: 
  teestream(ostream &o1, ostream &o2);
  
private :
  teebuf tbuf;   
}; 

