#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

#include <memory>
#include <string>

#include "axml.h"

using namespace std; 


class CommandLine
{
public: 
  CommandLine(int argc, char* argv[]);
  int getSeed() const {return seed; }
  string getConfigFileName() const {return configFileName; }
  string getAlnFileName() const{return alnFileName; }
  string getRunid() const {return runid; }
  string getTreeFile() const {return treeFile; }
  int getNumRunParallel() const {return runNumParallel; }

private: 
  int seed; 
  string configFileName; 
  string alnFileName; 
  string runid; 
  string treeFile; 
  string workDir;
  int runNumParallel; 

  void assertFileExists(string filename); 
  void parse(int argc, char *argv[]); 
  void printVersion();
  void printHelp();
}; 


#endif
