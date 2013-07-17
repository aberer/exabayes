#ifndef _CHECKPOINTABLE_H 
#define _CHECKPOINTABLE_H 

#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>

#define PRECISION (std::setprecision(std::numeric_limits<double>::digits10 + 2)) 


class Checkpointable
{
public: 
  Checkpointable()
    : DELIM('&')
    , checkpointIsBinary(false)
  {}

  virtual void readFromCheckpoint( std::ifstream &in )  = 0 ; 
  virtual void writeToCheckpoint( std::ofstream &out) = 0;   

  void getOfstream(std::string name, std::ofstream &result); 
  void getIfstream(std::string name, std::ifstream &result ); 

  template<typename T> void cWrite(std::ofstream &out, T toWrite); 
  template<typename T> T cRead(std::ifstream &in); 

private: 			// METHODS
  void readDelimiter(std::ifstream &in) ; 

private: 			// ATTRIBUTES 
  char DELIM; 
  bool checkpointIsBinary; 
}; 


template<typename T>
void Checkpointable::cWrite(std::ofstream &out, T toWrite)
{
  if(checkpointIsBinary)
    assert(0); 
  else 
    {
      out << std::scientific << PRECISION;
      out << toWrite << DELIM;         
    }
}


template<typename T>
T Checkpointable::cRead(std::ifstream &in )
{
  T result; 

  if(checkpointIsBinary)
    assert(0); 
  else 
    {
      in >> result; 
      // std::cout <<  "READ " << result << std::endl; 
      readDelimiter(in);
    }

  return result; 
}

#endif
