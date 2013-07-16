#ifndef _CHECKPOINTABLE_H 
#define _CHECKPOINTABLE_H 

#include <cassert>
#include <fstream>


class Checkpointable
{
public: 
  Checkpointable()
    : DELIM('&')
  { }

  virtual void readFromCheckpoint( std::ifstream &in )  = 0 ; 
  virtual void writeToCheckpoint( std::ofstream &out) const = 0;   
  
  void readDelimiter(std::ifstream &in) 
  {
    char c; 
    in >> c ; 
    assert(c == DELIM); 
  }
  
  char DELIM; 
  
}; 


#endif
