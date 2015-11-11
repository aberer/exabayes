#include "Checkpointable.hpp"

#include "GlobalVariables.hpp"


void Checkpointable::getOfstream(std::string name, std::ofstream &result)
{
  if(checkpointIsBinary)
    result.open(name, std::ios::binary); 
  else 
    result.open(name); 
}



void Checkpointable::getIfstream(std::string name, std::ifstream &result )
{
  if(checkpointIsBinary)
    result.open(name, std::ios::binary); 
  else 
    result.open(name); 
}



void Checkpointable::readDelimiter(std::ifstream &in)  
{
  char c; 
  in >> c ; 
  assert(c == DELIM); 
}


std::string Checkpointable::readString(std::ifstream &in )
{
  std::string result; 

  if(checkpointIsBinary)
    {
      int length = 0; 
      in.read((char*)&length, sizeof(int)); 
      
      result.resize(length) ;
      char *aString  = new char[length]; 
      in.read(aString, length * sizeof(char));       
      result = aString; 
      // std::cout << "read string >" << result << "< of length "<< result.size() << std::endl; 
    }
  else 
    {      
      getline(in, result, DELIM); 
    }

  return result; 
}

void Checkpointable::writeString(std::ofstream &out, std::string toWrite)
{
  if(checkpointIsBinary)
    {
      toWrite += '\0'; 
      int length = toWrite.size(); 
      out.write((char*)&length, sizeof(int)); 
      out.write(toWrite.c_str(), length * sizeof(char)); 
      // std::cout << "wrote string >"<< toWrite << "<  of length "  << toWrite.size() << std::endl; 
    }
  else  
    out << toWrite << DELIM; 
}
