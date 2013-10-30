#ifndef TREE_RESOURCE_HPP
#define TREE_RESOURCE_HPP

#include "InitializationResource.hpp"
class TreeAln; 

class TreeResource : public InitializationResource
{
public: 
  TreeResource(const TreeAln* tralnPtr);

  virtual std::vector<std::string> getTaxonNames(nat numTax) ; 
  virtual void fillAliasWgt(int *pos, nat length) ;   
  virtual std::tuple<int,int,double,int> getGlobalInfo() ;
  virtual void fillPartition(pInfo &partition, nat model)  ;
  virtual void fillAlnPart(unsigned char* ptr, nat length, nat &ctr) ; 
  virtual void fillParsVect(parsimonyNumber*& ptr, size_t &len, nat mult, nat model) ; 

  virtual bool isDataAreDistributed() { return true; }

private: 
  const TreeAln * _tralnPtr;
}; 

#endif
