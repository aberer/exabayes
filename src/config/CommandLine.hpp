#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

#include <memory>
#include <string>

#include "axml.h"

#include "Randomness.hpp"
#include "MemoryMode.hpp"


class CommandLine
{
public: 
  CommandLine(int argc, char* argv[]);
  randCtr_t getSeed() const; 
  std::string getConfigFileName() const {return configFileName; }
  std::string getAlnFileName() const{return alnFileName; }
  std::string getRunid() const {return runid; }
  std::string getTreeFile() const {return treeFile; }
  int getNumRunParallel() const {return runNumParallel; }
  std::string getWorkdir() const {return workDir; }
  void printVersion(bool toInfoFile);
  nat getNumChainsParallel() const {return chainNumParallel; }
  std::string getCheckpointId()const {return checkpointId; }
  bool isPerPartitionDataDistribution() const {return perPartitionDataDistribution; }
  void parseAlternative(int argc, char *argv[]); 

  MemoryMode getMemoryMode()const {return memoryMode ;  }

private: 			// METHODS

  void assertFileExists(std::string filename); 
  void parse(int argc, char *argv[]); 
  void printHelp();

private: 			// ATTRIBUTES
  randCtr_t seed; 
  std::string configFileName; 
  std::string alnFileName; 
  std::string runid; 
  std::string treeFile; 
  std::string workDir;
  int runNumParallel;   
  int chainNumParallel; 
  std::string checkpointId; 
  MemoryMode memoryMode; 
  bool perPartitionDataDistribution; 
  bool saveMemorySEV; 
}; 


#endif
