#ifndef BYTEFILERESOURCE
#define BYTEFILERESOURCE

#include "InitializationResource.hpp"

#include <fstream>

class ByteFileResource :  public InitializationResource
{
public: 
  ByteFileResource(std::string fileName);

  virtual void fillPartition(pInfo &partition, nat model) ; 
  virtual std::tuple<int,int,double,int> getGlobalInfo();
  virtual std::vector<std::string> getTaxonNames(nat numTax); 
  virtual void fillAliasWgt(int *pos, nat length); 
  virtual void fillAlnPart(unsigned char* ptr, nat length, nat &ctr);
  virtual void fillParsVect(parsimonyNumber*& ptr, size_t &len, nat mult, nat model); 

  virtual bool isDataAreDistributed() { return false; }
private: 
  std::ifstream byteFile; 
}; 

#endif
