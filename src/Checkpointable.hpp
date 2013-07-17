#ifndef _CHECKPOINTABLE_H 
#define _CHECKPOINTABLE_H 

#include <cassert>
#include <fstream>


class Checkpointable
{
public: 

  virtual void readFromCheckpoint( std::ifstream &in )  = 0 ; 
  virtual void writeToCheckpoint( std::ofstream &out) = 0;   
  
  void readDelimiter(std::ifstream &in) 
  {
    char c; 
    in >> c ; 
    assert(c == DELIM); 
  }
  
  static const char DELIM; 
  static const bool checkpointIsBinary; 

  static void getOfstream(std::string name, std::ofstream &result)
  {
    if(checkpointIsBinary)
      {
	result.open(name, std::ios::binary); 
	assert(0); 		// does not work that way 
      }
    else 
      {
	result.open(name); 
      }
  }

  static void getIfstream(std::string name, std::ifstream &result )
  {
    if(checkpointIsBinary)
      {
	result.open(name, std::ios::binary); 
	assert(0); 		// does not work that way 
      }
    else 
      result.open(name); 
  }
  
}; 


#endif
