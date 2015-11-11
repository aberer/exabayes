#ifndef  _CREDIBLE_SET_H
#define  _CREDIBLE_SET_H

#include "contrib/BipartitionExtractor.hpp" 


class CredibleSet
{
public: 
  CredibleSet(std::vector<std::string> file); 
  void printCredibleSet(std::string filename, double thresh)  ; 

private: 
  BipartitionExtractor _bipEx ; 
  // nat _totalTrees;   
}; 


#endif
