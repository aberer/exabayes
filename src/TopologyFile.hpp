#ifndef TOPOLOGY_FILE
#define TOPOLOGY_FILE

#include <sstream>
#include <string>
#include <iostream> 

#include "common.h"
#include "TreeAln.hpp"

class TopologyFile
{
public: 
  TopologyFile(std::string workdir, std::string runname, nat runid, nat couplingId); 
  void initialize(const TreeAln& traln, nat someId) const ;     
  void sample(const TreeAln &traln, nat gen) const ; 
  void finalize() const ; 
  void regenerate(std::string prevId, nat gen) ; 

private: 
  std::string fullFilename; 
  nat runid; 
  nat couplingId; 
}; 

#endif
