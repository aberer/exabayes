#ifndef _CHECKPOINTABLE_H 
#define _CHECKPOINTABLE_H 

#include <fstream>

class Checkpointable
{
public: 
  Checkpointable()
    : DELIM('&')
  { }

  virtual void readFromCheckpoint( std::ifstream &in )  = 0 ; 
  virtual void writeToCheckpoint( std::ofstream &out) const = 0;   
  
  char DELIM; 
  
}; 


#endif
