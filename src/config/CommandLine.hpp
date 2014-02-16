#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

#include <memory>
#include <string>

#include "axml.h"

#include "Randomness.hpp"
#include "MemoryMode.hpp"

#include "RunModes.hpp"


class CommandLine
{
public: 
  CommandLine();

  void initialize(  int argc, char **argv); 

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
  bool isSaveMemorySEV() const {return saveMemorySEV; }

  MemoryMode getMemoryMode()const {return memoryMode ;  }
  bool isDryRun() const {return dryRun; }
  RunModes getTreeInitRunMode() const ; 
  
  bool alnFileIsBinary() const; 

  std::string getSingleModel() const {return singleModel; }
  std::string getModelFile() const {return modelFile; }

  bool isQuiet() const {return quiet; }

private: 			// METHODS
  void assertFileExists(std::string filename); 
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
  bool dryRun; 
  std::string modelFile; 
  std::string singleModel; 
  bool quiet; 
}; 


#endif
