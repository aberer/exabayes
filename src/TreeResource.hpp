#ifndef TREE_RESOURCE_HPP
#define TREE_RESOURCE_HPP
#include <cassert>

#include "InitializationResource.hpp"
class TreeAln; 

class TreeResource : public InitializationResource
{
public: 
  TreeResource(const TreeAln* tralnPtr);
  virtual ~TreeResource(){}

  virtual std::vector<std::string> getTaxonNames(nat numTax) ; 
  virtual void fillAliasWgt(int *pos, nat length) ;   
  virtual std::tuple<int,int,double,int> getGlobalInfo() ;
  virtual void fillPartition(pInfo &partition, nat model)  ;
  virtual void fillAlnPart(unsigned char* ptr, nat length, nat &ctr) ; 

  virtual std::tuple<parsimonyNumber*,nat> fillParsVect( nat numTax, nat states, nat model); 

  virtual bool isDataAreDistributed() { return true; }
  
  virtual std::vector<double> getPartitionContributions(nat num) ; 
  
  virtual void initWeightsAndAln(TreeAln &traln)  ; 

private: 
  const TreeAln * _tralnPtr;
}; 

#endif
