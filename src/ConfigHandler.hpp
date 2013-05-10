#ifndef _CONFIGHANDLER_H
#define _CONFIGHANDLER_H

#include <string>

#include "ConfigReader.hpp"

using namespace std; 


class ConfigHandler
{
public: 
  ConfigHandler(string _configFileName) : configFileName(_configFileName) {};
  
  
  

private:   
  string configFileName; 

}; 


#endif
