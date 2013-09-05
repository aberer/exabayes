#include <iosfwd>

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "BipartitionExtractor.hpp" 


int main(int argc, char** argv)
{
  if(argc < 3) 
    {
      std::cout << "USAGE: ./extractBips id file[..]" << std::endl; 
      exit(-1); 
    }

  auto files = std::vector<std::string>();
  
  for(int i = 2; i < argc; ++i)
    {
      std::ifstream file(argv[i]);
      if(not file)
	{
	  std::cerr << "could not open file >" << argv[i] << "<" << std::endl; 
	  exit(-1); 
	}
      files.push_back(std::string(argv[i])); 
    }
  
  BipartitionExtractor  bipEx(files);

  std::string id = argv[1]; 

  bipEx.extractBipsNew();
  bipEx.printBipartitions(id);
  bipEx.printBipartitionStatistics(id); 
  bipEx.printFileNames(id);
  bipEx.printBranchLengths(id);

  return 0; 
}


