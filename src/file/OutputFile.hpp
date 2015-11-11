#ifndef OUTPUT_FILE_H 
#define OUTPUT_FILE_H 

#include <sys/stat.h>
#include <sstream>

class OutputFile
{
public: 
  static void rejectIfExists(std::string fileName); 
  void rejectIfNonExistant(std::string fileName); 
  // TODO portability 
  static std::string getFileBaseName(std::string workdir ); 
  std::string getFileName() const  {return fullFileName; }
  
  static bool directoryExists(std::string name); 

protected:   
  std::string fullFileName;   
}; 

#endif
