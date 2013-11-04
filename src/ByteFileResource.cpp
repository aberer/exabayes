#include <cassert>
#include <cstring>

#include "ByteFileResource.hpp"
#include "GlobalVariables.hpp"


ByteFileResource::ByteFileResource(std::string fileName)
  : byteFile{fileName, std::ios::binary}
{
  // parse the "magic number"
  char tmp[7]; 
  memset(&tmp, 0 , sizeof(char) * 7 ); 
  byteFile.read(tmp, 6 * sizeof(char)); 
  if( std::string(tmp).compare("BINARY") != 0 )
    {
      tout << "expected different file start (maybe reparse the binary file). Got >"  << tmp << "<"<< std::endl; 
      assert(0); 
    }
}


template<typename T>
static void byteRead(std::ifstream& in, T* result, nat num)
{
  in.read((char*)result, sizeof(T) * num ); 
}



std::tuple<int,int,double,int> ByteFileResource::getGlobalInfo()
{
  auto mxTips = int{0} ;
  auto numberOfPartitions = int{0}; 
  auto gappyness = double{0.};
  auto origLength = int{0}; 

  byteRead(byteFile, &mxTips,1); // 1
  byteRead(byteFile, &numberOfPartitions, 1);	      // 2
  byteRead(byteFile, &gappyness,1);		      // 4
  byteRead(byteFile, &origLength, 1); // 3
  
  
  return std::make_tuple(mxTips, numberOfPartitions, gappyness, origLength);
}


void ByteFileResource::fillAliasWgt(int *pos, nat length)
{
  byteRead(byteFile, pos, length);
}

std::vector<std::string> ByteFileResource::getTaxonNames(nat numTax)
{  
  auto nameMap = std::vector<std::string>{}; 
  for(nat i = 1; i <= numTax; i++)
    {
      int 
  	len = 0;

      byteRead(byteFile, &len, 1); 
      auto tmp = new char[len]; 
      byteRead(byteFile, tmp , len); 
      nameMap.push_back(tmp);
      delete [] tmp ; 
    }  
  return nameMap; 
}


void ByteFileResource::fillPartition(pInfo &partition, nat model)
{
  byteRead(byteFile, &partition.states, 1); 
  byteRead(byteFile, &partition.maxTipStates, 1); 
  byteRead(byteFile, &partition.lower, 1); 
  byteRead(byteFile, &partition.upper, 1); 
  byteRead(byteFile, &partition.width, 1); 
  byteRead(byteFile, &partition.dataType, 1); 
  byteRead(byteFile, &partition.protModels, 1); 
  byteRead(byteFile, &partition.protFreqs, 1); 
  byteRead(byteFile, &partition.nonGTR, 1); 
      
  int len = 0; 
  byteRead(byteFile, &len, 1); 
  partition.partitionName = (char*) exa_calloc( len, sizeof(char));
  byteRead(byteFile, partition.partitionName, len); 

}


void ByteFileResource::fillAlnPart(unsigned char* ptr, nat length, nat &ctr)
{
  byteRead(byteFile, ptr, length); 
  ++ctr; 
}


void ByteFileResource::fillParsVect(parsimonyNumber*& ptr, size_t &len, nat mult, nat model)
{
  byteRead(byteFile,&len , 1); 
  nat numBytes =  mult * len; 
  ptr = (parsimonyNumber*)exa_malloc_aligned( numBytes * sizeof(parsimonyNumber));
  memset(ptr, 0 , sizeof(parsimonyNumber) * numBytes) ;
  byteRead(byteFile, ptr, numBytes); 
}


std::vector<double> ByteFileResource::getPartitionContributions(nat num) 
{
  auto result = std::vector<double>{}; 
  for(nat i = 0; i < num; ++i)
    {
      double elem = 0; 
      byteRead(byteFile, &elem,1); 
      result.push_back(elem); 
    }
  return result; 
} 
