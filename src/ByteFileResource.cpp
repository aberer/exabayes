#include <cassert>
#include <cstring>

#include "ByteFileResource.hpp"
#include "GlobalVariables.hpp"


ByteFileResource::ByteFileResource(std::string fileName, ParallelSetup  pl)
  : _byteFile{fileName, std::ios::binary}
  , _pl(pl)
{
  auto tmp = std::vector<char> (7, '\0'); 
  _byteFile.read(tmp.data(), 6 * sizeof(char)); 
  auto asString = std::string(begin(tmp), end(tmp));
  if( asString.compare("BINARY") == 0 )
    {
      tout << "expected different file start (maybe reparse the binary file). Got >"  << asString << "<"<< std::endl; 
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

  byteRead(_byteFile, &mxTips,1); // 1
  byteRead(_byteFile, &numberOfPartitions, 1);	      // 2
  byteRead(_byteFile, &gappyness,1);		      // 4
  byteRead(_byteFile, &origLength, 1); // 3

  return std::make_tuple(mxTips, numberOfPartitions, gappyness, origLength);
}


void ByteFileResource::fillAliasWgt(int *pos, nat length)
{
  byteRead(_byteFile, pos, length);
}

std::vector<std::string> ByteFileResource::getTaxonNames(nat numTax)
{  
  auto nameMap = std::vector<std::string>{}; 
  for(nat i = 1; i <= numTax; i++)
    {
      int 
  	len = 0;

      byteRead(_byteFile, &len, 1); 
      auto tmp = new char[len]; 
      byteRead(_byteFile, tmp , len); 
      nameMap.push_back(tmp);
      delete [] tmp ; 
    }  
  return nameMap; 
}


void ByteFileResource::fillPartition(pInfo &partition, nat model)
{
  byteRead(_byteFile, &partition.states, 1); 
  byteRead(_byteFile, &partition.maxTipStates, 1); 
  byteRead(_byteFile, &partition.lower, 1); 
  byteRead(_byteFile, &partition.upper, 1); 
  byteRead(_byteFile, &partition.width, 1); 
  byteRead(_byteFile, &partition.dataType, 1); 
  byteRead(_byteFile, &partition.protModels, 1); 
  byteRead(_byteFile, &partition.protFreqs, 1); 
  byteRead(_byteFile, &partition.nonGTR, 1); 
      
  int len = 0; 
  byteRead(_byteFile, &len, 1); 
  partition.partitionName = (char*) exa_calloc( len, sizeof(char));
  byteRead(_byteFile, partition.partitionName, len); 

}


void ByteFileResource::fillAlnPart(unsigned char* ptr, nat length, nat &ctr)
{
  byteRead(_byteFile, ptr, length); 
  ++ctr; 
}



static void printBits(parsimonyNumber num )
{
  for(nat i = 0; i < 32; ++i)
    {
      if( ( num  & ( 1 << i ) )  != 0  ) 
	tout << "1" ; 
      else 
	tout << "0";
    }
}


std::tuple<parsimonyNumber*,nat> ByteFileResource::fillParsVect( nat numNode, nat states, nat model)
{
  parsimonyNumber* ptr; 
  nat len = 0; 
  nat mult = numNode * states; 
 
#if HAVE_PLL != 0
  // old stuff, where nothing is distributed 
  byteRead(_byteFile,&len , 1); 

  assert(len != 0 ); 

  nat numBytes =  mult * len; 
  ptr = (parsimonyNumber*)exa_malloc_aligned( numBytes * sizeof(parsimonyNumber));
  // memset(ptr, 0 , sizeof(parsimonyNumber) * numBytes) ;

  // tout << "reading  " <<  numBytes << std::endl; 
  byteRead(_byteFile, ptr, numBytes); 
#else 
  // data distribution is not necessary in line with partition
  // distribution scheme. But this will not hurt, unless we are trying
  // to deal with ~ 1e4 partitions

  nat totalLen = 0; 
  byteRead(_byteFile,&totalLen , 1); 

  assert(totalLen != 0 ); 

  nat numBytes =  mult * totalLen; 

  auto allArray = std::vector<parsimonyNumber> (numBytes, 0); 
  byteRead(_byteFile, allArray.data(), numBytes); 
  // tout << "model=" << model << "\tlength=" << totalLen << "\treading  " <<  numBytes << std::endl; 

  nat rank = _pl.getRankInChainBatch(); 
  nat size = _pl.getProcessesPerChainBatch();

  len = (totalLen / size); 
  if( rank < totalLen % size )
    ++len; 
  auto lenUnpadded = len; 
  
  // tout << "totallen=" << totalLen << std::endl; 

  // could be avx-specific, but it's hardly worth the coding overhead... 
  const int intsPerType = 8; 

  // correct the padding 
  if( len % intsPerType != 0 ) 
    len +=  intsPerType - ( len % intsPerType ) ; 

  ptr = (parsimonyNumber*) exa_malloc_aligned( ( len * mult )  * sizeof(parsimonyNumber)); 
  memset(ptr, std::numeric_limits<parsimonyNumber>::max(),( len * mult )  * sizeof(parsimonyNumber)) ; 

  // redistribute the data 
  nat numTax = numNode / 2 ;
  assert(numNode % 2  == 0); 	// meh

  for(nat i = 0; i <  2 * numTax   ; ++ i )
    {
      // tout << "taxon "<< i << std::endl; 
      for(nat k = 0; k < states ; ++k )
	{
	  // tout << "state "<< k <<  " " ; 
	  auto curStart = &(allArray.at( totalLen * states * (i )  + ( totalLen * k  )  ))  ; 
	  auto localPtr = ptr + len * states * (i ) + (len * k ); 
	    
	  nat ctr = 0; 
	  for( nat j = 0 ; j < totalLen; ++j )
	    {
	      if(j % size == rank)
		{
		  // printBits(curStart[j]) ; 
		  // tout << ","; 
		  localPtr[ctr++] = curStart[j]; 
		}
	    }
	  // tout << std::endl; 
	  assert(ctr == lenUnpadded); 
	}
      // tout << std::endl; 
    }
#endif

  // tout << "myLen=" << len << std::endl; 

  return std::make_tuple(ptr, len);
}


std::vector<double> ByteFileResource::getPartitionContributions(nat num) 
{
  auto result = std::vector<double>{}; 
  for(nat i = 0; i < num; ++i)
    {
      double elem = 0; 
      byteRead(_byteFile, &elem,1); 
      result.push_back(elem); 
    }
  return result; 
} 
