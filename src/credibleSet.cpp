#include "CredibleSet.hpp"

#include <sstream>
#include <iosfwd>

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

int main(int argc, char **argv)
{
  if(argc != 4)
    {
      std::cerr << "./credibleSet id  credibleInterval treeFile" <<  std::endl; 
      std::cerr << "\twhere credibleInterval is the percentile (e.g., 50 for the median)\n"
		<< "\tof trees that are the most frequent" << std::endl;
      exit(-1);       
    }

  auto id = std::string(argv[1]); 
  auto cIString = std::string(argv[2]) ; 

  std::istringstream  os{cIString}; 
  double ci = 0.; 
  os >> ci; 
  assert(ci > 1); 
  ci /= 100.; 

  auto file = std::string(argv[3]) ; 
  auto cs = CredibleSet(file); 

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_credibleSet." << id  ; 

  if(std::ifstream(ss.str()))
    {
      std::cerr << "The file " << ss.str() << " already exists (possibly left over from a previous run). Please\n"
		<< "choose a different run-id." << std::endl; 
      exit(-1); 
    }

  cs.printCredibleSet(ss.str(), ci); 

  return 0; 
}
