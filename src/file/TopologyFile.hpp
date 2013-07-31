#ifndef TOPOLOGY_FILE
#define TOPOLOGY_FILE

#include <sstream>
#include <string>
#include <iostream> 

#include "common.h"
#include "TreeAln.hpp"
#include "OutputFile.hpp"

class TopologyFile : public OutputFile
{
public: 
  TopologyFile(std::string workdir, std::string runname, nat runid, nat couplingId); 
  void initialize(const TreeAln& traln, nat someId)  ;     
  void sample(const TreeAln &traln, nat gen)  ; 
  void finalize()  ; 
  void regenerate(std::string workdir, std::string prevId, nat gen) ; 
  void verifyNonExistance();

private: 
  nat runid; 
  nat couplingId; 
}; 

#endif
