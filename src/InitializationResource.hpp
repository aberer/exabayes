#ifndef INITRESOURCE_HPP
#define INITRESOURCE_HPP

#include "common.h"
#include "axml.h"

#include "TreeAln.hpp"

#include <algorithm>

class InitializationResource 
{
public: 
  virtual ~InitializationResource(){}
  virtual std::vector<std::string> getTaxonNames(nat numTax) = 0; 
  virtual void fillAliasWgt(int *pos, nat length) = 0;   
  virtual std::tuple<int,int,double,int> getGlobalInfo() = 0;
  virtual void fillPartition(pInfo &partition, nat model)  = 0;
  virtual void fillAlnPart(unsigned char* ptr, nat length, nat &ctr) = 0; 
  virtual std::tuple<parsimonyNumber*,nat> fillParsVect( nat numTax, nat states, nat model) = 0; 
  virtual std::vector<double> getPartitionContributions(nat num)  = 0; 
  virtual bool isDataAreDistributed() = 0; 
  virtual void initWeightsAndAln(TreeAln &traln)  = 0; 
}; 



#endif
