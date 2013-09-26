#include "Serializable.hpp"
#include "GlobalVariables.hpp"


void Serializable::getOfstream(std::string name, std::ofstream &result)
{
  if(checkpointIsBinary)
    result.open(name, std::ios::binary); 
  else 
    result.open(name); 
}



void Serializable::getIfstream(std::string name, std::ifstream &result )
{
  if(checkpointIsBinary)
    result.open(name, std::ios::binary); 
  else 
    result.open(name); 
}



void Serializable::readDelimiter(std::istream &in)  
{
  char c; 
  in >> c ; 
  assert(c == DELIM); 
}


std::string Serializable::readString(std::istream &in )
{
  auto result = std::string{}; 

  if(checkpointIsBinary)
    {
      int length = 0; 
      in.read((char*)&length, sizeof(int)); 
      
      result.resize(length) ;
      char *aString  = new char[length]; 
      in.read(aString, length * sizeof(char));       
      result = aString; 
    }
  else 
    {      
      getline(in, result, DELIM); 
    }
  return result; 
}

void Serializable::writeString(std::ostream &out, std::string toWrite) const 
{  
  if(checkpointIsBinary)
    {
      toWrite += '\0'; 
      int length = toWrite.size(); 
      out.write((char*)&length, sizeof(int)); 
      out.write(toWrite.c_str(), length * sizeof(char)); 
    }
  else  
    {
      out << toWrite << DELIM; 
    }
}
