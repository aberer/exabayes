#ifndef OUTPUT_FILE_H 
#define OUTPUT_FILE_H 


#include "ParallelSetup.hpp"

class OutputFile
{
public: 
  void rejectIfExists(std::string fileName)
  {
    if( std::ifstream(fileName) ) 
      {
	std::cerr << std::endl <<  "File " << fileName << " already exists (probably \n"
		  << "from previous run). Please choose a new run-id or remove previous output files. " << std::endl; 
	ParallelSetup::genericExit(-1);
    }
  }

  
  void rejectIfNonExistant(std::string fileName)
  {
    if( not std::ifstream(fileName) )      
      {
	std::cerr << "Error: could not find file from previous run. \n"
		  << "The assumed name of this file was >" <<  fileName << "<. Aborting." << std::endl; 
	ParallelSetup::genericExit(0); 
      }
  }


protected:   
  std::string fullFileName;   
}; 

#endif
