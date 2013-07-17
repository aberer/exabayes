#include "Checkpointable.hpp"

#include "GlobalVariables.hpp"


void Checkpointable::getOfstream(std::string name, std::ofstream &result)
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



void Checkpointable::getIfstream(std::string name, std::ifstream &result )
{
  if(checkpointIsBinary)
    {
      result.open(name, std::ios::binary); 
      assert(0); 		// does not work that way 
    }
  else 
    result.open(name); 
}



void Checkpointable::readDelimiter(std::ifstream &in)  
{
  char c; 
  in >> c ; 
  assert(c == DELIM); 
}


template<>
void Checkpointable::cWrite<std::string>(std::ofstream &out, std::string toWrite)
{
  if(checkpointIsBinary)
    {
      assert(0); 
    }
  else  
    {
      out << toWrite << DELIM; 
    }
}



template<>
std::string Checkpointable::cRead(std::ifstream &in )
{
  std::string result; 

  if(checkpointIsBinary)
    {
      assert(0); 
    }
  else 
    {
      getline(in, result, DELIM); 
    }

  return result; 
}
