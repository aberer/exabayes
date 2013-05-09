#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

#include <memory>

#include "axml.h"
// #include "bayes.h"


class CommandLine
{
public: 
  CommandLine(int argc, char* argv[]);  

  analdef* getAdef() const {return adef.get(); } 
  int getSeed() const {return seed; }

private: 
  int seed; 
  unique_ptr<analdef>  adef; 
  void parse(int argc, char *argv[]); 
  void initAdef();
  void printVersion();
  void printHelp();
}; 


#endif
