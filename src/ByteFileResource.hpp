#ifndef BYTEFILERESOURCE
#define BYTEFILERESOURCE


#include "ParallelSetup.hpp"
#include "InitializationResource.hpp"

#include <memory>
#include <fstream>

class ByteFileResource :  public InitializationResource
{
public: 
  ByteFileResource(std::string fileName, std::shared_ptr<ParallelSetup> pl );
  virtual ~ByteFileResource(){}

  virtual void fillPartition(pInfo &partition, nat model) ; 
  virtual std::tuple<int,int,double,int> getGlobalInfo();
  virtual std::vector<std::string> getTaxonNames(nat numTax); 
  virtual void fillAliasWgt(int *pos, nat length); 
  virtual void fillAlnPart(unsigned char* ptr, nat length, nat &ctr);
  virtual std::tuple<parsimonyNumber*,nat> fillParsVect(nat numTax, nat states, nat model); 

  virtual std::vector<double> getPartitionContributions(nat num) ; 
  virtual void initWeightsAndAln(TreeAln &traln) {assert(0) ; }
  virtual bool isDataAreDistributed() { return false; } 
private: 
  std::ifstream _byteFile; 
  std::shared_ptr<ParallelSetup> _pl; 		// non-owning
}; 

#endif
