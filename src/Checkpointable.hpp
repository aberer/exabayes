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
    , checkpointIsBinary(true)
  {}

  virtual void readFromCheckpoint( std::ifstream &in )  = 0 ; 
  virtual void writeToCheckpoint( std::ofstream &out) = 0;   

  void getOfstream(std::string name, std::ofstream &result); 
  void getIfstream(std::string name, std::ifstream &result ); 

  template<typename T> void cWrite(std::ofstream &out, T &toWrite); 
  template<typename T> T cRead(std::ifstream &in); 

  std::string readString(std::ifstream &in ); 
  void writeString(std::ofstream &out, std::string toWrite); 

private: 			// METHODS
  void readDelimiter(std::ifstream &in) ; 

private: 			// ATTRIBUTES 
  char DELIM; 
  bool checkpointIsBinary; 
}; 


template<typename T>
void Checkpointable::cWrite(std::ofstream &out, T& toWrite)
{
  static_assert(not std::is_same<T, std::string>::value, "Do NOT use the cWrite funciton with strings (there is a specific function for that)");

  if(checkpointIsBinary)
    {      
      out.write((char*)&toWrite, sizeof(T)); 
      // std::cout << "WROTE " << toWrite << std::endl;
    }
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
