#include <cassert>
#include <cstring>

#include <numeric>
#include <algorithm>

#include "ByteFileResource.hpp"
#include "GlobalVariables.hpp"

#include "TreeAln.hpp"

// #define DEBUG_MSG


ByteFileResource::ByteFileResource(std::string fileName, std::shared_ptr<ParallelSetup> pl)
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
static void byteRead(std::ifstream& in, T* result, size_t num)
{
  assert(num <= num * sizeof(T) ); 
  
  in.read((char*)result, sizeof(T) * num ); 
  auto bytesRead = in.gcount();
  assert(bytesRead == std::streamsize(sizeof(T) * num));
}


std::tuple<int,int,double,int> ByteFileResource::getGlobalInfo()
{
#ifdef DEBUG_MSG
  tout << "getting global info..  "; 
#endif
  auto mxTips = int{0} ;
  auto numberOfPartitions = int{0}; 
  auto gappyness = double{0.};
  auto origLength = int{0}; 

  byteRead(_byteFile, &mxTips,1); // 1
  byteRead(_byteFile, &numberOfPartitions, 1);	      // 2
  byteRead(_byteFile, &gappyness,1);		      // 4
  byteRead(_byteFile, &origLength, 1); // 3

#ifdef DEBUG_MSG
  tout << "done  " << std::endl; 
#endif

  return std::make_tuple(mxTips, numberOfPartitions, gappyness, origLength);
}

#if HAVE_PLL == 0
void ByteFileResource::fillAliasWgt(TreeAln& traln)
{
#ifdef DEBUG_MSG
  tout << "filling alias wgt.." ; 
#endif

  // reset the stream to where the weights are located 
  auto posToRestore =  _byteFile.tellg();
  _byteFile.seekg(_weightPos, ios_base::beg);

  // now distribute the weights 
  auto &tr = traln.getTrHandle(); 

  // all processes must know, how many weigths they receive and what
  // displacement they have in the array
  auto countsPerProc = std::vector<int>(processes);
  for(nat model = 0; model < traln.getNumberOfPartitions() ;++model)
    {
      auto &partition = traln.getPartition(model);
      auto length = partition.upper - partition.lower; 
      if(tr.manyPartitions)
	{
	  auto currentRank = tr.partitionAssignment[model]; 
	  countsPerProc[currentRank] += length; 
	}
      else 
	{
	  for(int j = partition.lower; j < partition.upper ; ++j)
	    ++countsPerProc[j % processes]; 
	}
    }

  int sumOfSites = std::accumulate(begin(countsPerProc), end(countsPerProc), 0 ); 
  assert(sumOfSites == tr.originalCrunchedLength); 

  auto displ = std::vector<int>{}; 
  displ.push_back(0);
  for(int i = 1 ; i < processes; ++i)
    displ.push_back(displ[i-1] + countsPerProc[i-1]);
  auto myData = std::vector<int>(countsPerProc[processID]); 

  // BAD: we cannot use the parallel setup at this point, because we
  // do not have access to it in the TreeAln operator. 
  if(processID == 0)
    {
      // read the data once 
      auto allData = std::vector<int>(tr.originalCrunchedLength);
      byteRead(_byteFile, allData.data(), tr.originalCrunchedLength);

      // reorder the data according to the scheme 
      auto allDataReordered = std::vector<int>(tr.originalCrunchedLength); 
      auto indexPerProc = displ; 
      for(nat model = 0; model < traln.getNumberOfPartitions() ; ++model)
	{			
	  auto &partition = traln.getPartition(model);
	  
	  // partition distribution 
	  if(tr.manyPartitions)
	    {
	      auto currentRank = tr.partitionAssignment[model]; 
	      std::copy(begin(allData) + partition.lower, begin(allData) + partition.upper, begin(allDataReordered) + indexPerProc[currentRank]);
	      indexPerProc[currentRank] += partition.upper - partition.lower; 
	    }
	  // cyclic distribution 
	  else 
	    {
	      for(int i = partition.lower; i < partition.upper; ++i)
		{
		  auto currentRank = i % processes; 
		  allDataReordered[indexPerProc[currentRank]] = allData[i]; 
		  ++indexPerProc[currentRank]; 
		}
	    }
	}

      // the barrier and all that shrinking is a bit excessive. I fear
      // it's necessary for whole genome files
      
      allData.clear(); 
      allData.shrink_to_fit();

      // inform your peer processes 
      MPI_Scatterv(allDataReordered.data(), countsPerProc.data(), displ.data(), MPI_INT, myData.data(), countsPerProc[processID], MPI_INT, 0, comm);
      allDataReordered.clear(); 
      allDataReordered.shrink_to_fit();
      MPI_Barrier(comm);
    }
  else 
    {
      MPI_Scatterv(NULL, countsPerProc.data(), displ.data(), MPI_INT, myData.data(), countsPerProc[processID], MPI_INT, 0, comm);
      MPI_Barrier(comm);
    }

  // now every rank has exactly their data. We just have to sort it
  // back in
  auto iter = begin(myData); 
  for(nat model = 0; model < traln.getNumberOfPartitions(); ++model)
    {
      auto &partition = traln.getPartition(model);
      if(partition.width > 0 )
	{
	  std::copy(iter, iter + partition.width, partition.wgt); 
	  iter += partition.width; 
	}
    }
  assert(iter == end(myData));

  // reset the stream to the position, we've previously been 
  _byteFile.seekg(posToRestore,ios_base::beg); 

#ifdef DEBUG_MSG
  tout << "done" << std::endl; 
#endif
}
#else 
void ByteFileResource::fillAliasWgt(TreeAln& traln)
{
  // reset the stream to where the weights are located 
  auto posToRestore =  _byteFile.tellg();
  _byteFile.seekg(_weightPos, ios_base::beg);

  auto &tr = traln.getTrHandle(); 
  auto length = tr.originalCrunchedLength; 
  auto pos = new int[length]; 
  byteRead(_byteFile, pos, length);

  for(nat model = 0 ; model < traln.getNumberOfPartitions() ; ++model)
    {
      auto &partition = traln.getPartition(model); 
      // partition.wgt[]
      nat ctr = 0; 
      for(int i = partition.lower; i < partition.upper; ++i)
	partition.wgt[ctr++] = pos[ i ]; 
      // std::copy(pos + partition.lower, pos + partition.upper ,  partition.wgt ); 
    }
  
  // reset the stream to the position, we've previously been 
  _byteFile.seekg(posToRestore,ios_base::beg); 
}
#endif


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
#ifdef DEBUG_MSG
  tout << "filling partition " << model << "..."; 
#endif

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

#ifdef DEBUG_MSG
  tout << "done" << std::endl; 
#endif
}


void ByteFileResource::fillAlnPart(unsigned char* ptr, nat length, nat &ctr)
{
#ifdef DEBUG_MSG
  tout << "filling aln part " << ctr  <<  "..."; 
#endif

  byteRead(_byteFile, ptr, length); 
  ++ctr; 

#ifdef DEBUG_MSG
  tout << "done"<< std::endl; 
#endif
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
#ifdef DEBUG_MSG
  tout << "filling pars vect for model " << model << "..."; 
#endif

  parsimonyNumber* ptr; 
  nat len = 0; 
  size_t mult = numNode * states; 
 
#if HAVE_PLL != 0
  // old stuff, where nothing is distributed 
  byteRead(_byteFile,&len , 1); 

  assert(len != 0 ); 

  size_t numBytes =  mult * len; 
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

  // nat rank = _pl.getRankInChainBatch(); 
  nat rank = _pl->getRankInData();
  nat size = _pl->getSizeInData(); 

  // nat size = _pl->.getProcessesPerChainBatch();
  // nat size = _pl->.get

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
  assert(numNode % 2  == 0);    // meh

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

#ifdef DEBUG_MSG
  tout << "done"  << std::endl; 
#endif

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


void ByteFileResource::markPosAndSkipWeights(int len)
{
  _weightPos = _byteFile.tellg();
  _byteFile.seekg(len * sizeof(int), ios_base::cur);
} 

