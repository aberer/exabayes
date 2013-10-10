#include <sstream>
#include <iosfwd>

int NUM_BRANCHES; 	

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "ConsensusTree.hpp"


int main(int argc, char **argv)
{
  // BAD
  NUM_BRANCHES = 1; 

  if(argc != 4 ) 
    {
      std::cout << "Usage: ./consense id threshold file" << std::endl; 
      std::cout << "\twhere threshold is a consensus threshold (between 50 and 100 or MRE\n"
		<< "for a greedily refined majority-rule consensus tree)" << std::endl; 
      exit(-1); 
    }
  
  auto threshold = 0.; 
  auto thresholdString = std::string(argv[2]); 

  std::istringstream helper(thresholdString); 
  bool isMre = false; 
  if(thresholdString.compare("MRE") == 0 || thresholdString.compare("mre") == 0)
    {
      threshold = 50; 
      isMre = true; 
    }
  else 
    helper >> threshold; 
  
  auto id = std::string{argv[1]}; 
  auto file = argv[3]; 
  
  assert(threshold > 1); 
  threshold /= 100.; 

  auto ct = ConsensusTree(file); 
  auto result = ct.getConsensusTreeString(threshold,isMre);

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_consensusTree." << id; 

  if(std::ifstream(ss.str()))
    {
      std::cerr << std::endl << "File " << ss.str() << " already exists (probably \n"
		<< "left over from a previous run). Please choose a new run-id or remove\n"
		<< "previous output files." << std::endl; 
      exit(-1); 
    }

  std::ofstream outfile(ss.str()); 
  outfile << result << std::endl; 

  std::cout << "Printed consensus tree to " << ss.str() << std::endl; 

  return 0; 
}
