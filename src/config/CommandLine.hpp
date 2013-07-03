#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

#include <memory>
#include <string>

#include "axml.h"

#include "Randomness.hpp"

using namespace std; 


class CommandLine
{
public: 
  CommandLine(int argc, char* argv[]);
  randCtr_t getSeed() const; 
  string getConfigFileName() const {return configFileName; }
  string getAlnFileName() const{return alnFileName; }
  string getRunid() const {return runid; }
  string getTreeFile() const {return treeFile; }
  int getNumRunParallel() const {return runNumParallel; }
  string getWorkdir() const {return workDir; }
  void printVersion(bool toInfoFile);
  // int getChainNumParallel() const { return chainNumParallel; }

  void parseAlternative(int argc, char *argv[]); 

private: 
  // int seed; 
  randCtr_t seed; 
  string configFileName; 
  string alnFileName; 
  string runid; 
  string treeFile; 
  string workDir;
  int runNumParallel; 
  int chainNumParallel; 

  void assertFileExists(string filename); 
  void parse(int argc, char *argv[]); 
  void printHelp();
}; 


#endif
