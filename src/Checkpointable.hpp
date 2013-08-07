#ifndef _CHECKPOINTABLE_H 
#define _CHECKPOINTABLE_H 

#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>
#include "common.h"


class Checkpointable
{
public: 
  Checkpointable()
    : DELIM('&')
    , checkpointIsBinary(true)
  {}

  virtual void readFromCheckpoint( std::istream &in )  = 0 ; 
  virtual void writeToCheckpoint( std::ostream &out) const = 0;   

  void getOfstream(std::string name, std::ofstream &result); 
  void getIfstream(std::string name, std::ifstream &result ); 

  template<typename T> void cWrite(std::ostream &out, const T &toWrite) const; 
  template<typename T> T cRead(std::istream &in); 

  std::string readString(std::istream &in ); 
  void writeString(std::ostream &out, std::string toWrite) const; 

private: 			// METHODS
  void readDelimiter(std::istream &in) ; 

private: 			// ATTRIBUTES 
  char DELIM; 
  bool checkpointIsBinary; 
}; 


template<typename T>
void Checkpointable::cWrite(std::ostream &out, const T& toWrite) const 
{
  static_assert(not std::is_same<T, std::string>::value, "Do NOT use the cWrite funciton with strings (there is a specific function for that)");

  if(checkpointIsBinary)
    {      
      out.write((char*)&toWrite, sizeof(T)); 
      // std::cout << "WROTE " << toWrite << std::endl;
    }
  else 
    {
      out << std::scientific << MAX_SCI_PRECISION; 
      out << toWrite << DELIM; 
    }
}


template<typename T>
T Checkpointable::cRead(std::istream &in )
{
  T result; 

  static_assert(not std::is_same<T, std::string>::value, "Do NOT use the cRead funciton with strings (there is a specific function for that)");

  if(checkpointIsBinary)
    {
      in.read((char*)&result, sizeof(T)); 
      // std::cout << "READ " << result << std::endl; 
    }
  else 
    {
      in >> result; 
      // std::cout <<  "READ " << result << std::endl; 
      readDelimiter(in);
    }

  return result; 
}

#endif
