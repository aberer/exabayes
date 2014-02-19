#include <cassert>
#include "GlobalVariables.hpp"


template<typename T> 
auto  LocalComm::scatterVariableKnownLength( std::vector<T> allData, std::vector<int> &countsPerProc, std::vector<int> &displPerProc, int root)   
  -> std::vector<T>
{
  // std::cout << SyncOut() << "scatterVariableKnownLength" << SHOW(root) << std::endl; 
  int myRank = getRank(); 
  int myCol = getColor();
  int rootIdx = getIdx(myCol, root); 
  int myIdx = getIdx(); 
  auto &rootMsg = _messages.at(rootIdx); 
  auto &myMsg = _messages.at(myIdx); 
  
  auto result = std::vector<T>(countsPerProc[myRank],0); 
  
  if(myRank == root)
    rootMsg.post(allData);
  else 
    myMsg.incrementGeneration(true); 

  auto myGen = myMsg.getGeneration();

  // everybody reads 
  while( rootMsg.getGeneration() != myGen) 
    ; 
  
  rootMsg.readMessage( result.data(), displPerProc[myRank], countsPerProc[myRank]); 
  
  while(myRank == root && not rootMsg.clearMessage(size()))
    ;

  return result; 
} 


template<typename T>
auto LocalComm::gatherVariableKnownLength(std::vector<T> myData, std::vector<int> &countsPerProc, std::vector<int> &displPerProc , int root)   
  -> std::vector<T>
{
  // TODO, check for zero-messages and skip 
  auto myRank = getRank(); 
  auto myCol = getColor(); 
  auto myIdx = getIdx(); 

  auto result = std::vector<T>( std::accumulate(begin( countsPerProc), end( countsPerProc), 0),
				0); 
  auto resultIter = begin(result); 

  auto &msg = _messages.at(myIdx); 
  
  assert(int(myData.size()) == countsPerProc[myRank]) ;
  msg.post(myData);
  auto myGen =  msg.getGeneration();

  if(myRank == root)
    {
      nat gatheredTotal = 0; 
      for(int i = 0; i < size(); ++i)
	{
	  auto hisIdx = getIdx(myCol, i); 

	  auto &hisMsg = _messages.at(hisIdx); 

	  while( hisMsg.getGeneration() != myGen )
	    ; 

	  auto msgSiz = hisMsg.readMessage( &(*resultIter), 0, countsPerProc[i]); 
	  resultIter += msgSiz; 
	  gatheredTotal += msgSiz; 
	}

      assert(gatheredTotal == result.size()); 

      auto wasCleared = msg.clearMessage(1); 
      assert(wasCleared); 
    }
  else 
    {
      while(not msg.clearMessage(1)) // only root reads message 
	; 
    }

  return result; 
} 


template<typename T> 
auto LocalComm::broadcast(std::vector<T> array, int root ) 
  -> std::vector<T>
{
  int myRank = getRank();
  int myIdx = getIdx();
  int myCol = getColor();
  int rootIdx = getIdx(myCol, root); 

  auto& rootMsg = _messages.at(rootIdx);   
  auto &msg = _messages.at(myIdx); 

  if(getRank() == root)
    rootMsg.post(array);
  else 
    msg.incrementGeneration(true);

  auto myGen = msg.getGeneration();
    
  if(getRank() == root)
    {
      auto num = size() -1  ; 
      while(not msg.clearMessage(num ))
	; 
    }
  else 
    {
      while(myGen != rootMsg.getGeneration())
	; 
      rootMsg.readMessage(array.data(), 0, array.size());
    }

  return array; 
}



template<typename T>
std::vector<T> LocalComm::reduce(std::vector<T> data, int root )
{
  int myIdx = getIdx();
  int myCol = getColor();
  int myRank = getRank(); 
  assert(root < size());

  auto &msg = _messages.at(myIdx) ; 

  assert(root < size()); 
  assert(myRank < size() ); 
  
  msg.post(data); 
  auto myGen = msg.getGeneration(); 
  
  auto result = std::vector<T>(data.size(),0);
  
  if(myRank == root)
    {
      for(int i = 0; i < size() ;++i)
	{
	  int hisIdx = getIdx(myCol, i); 
	  auto &hisMsg = _messages.at(hisIdx); 
	  while( myGen !=  hisMsg.getGeneration())
	    ; 
	  auto hisData = hisMsg.getMessage<T>();

	  assert(hisData.size() == result.size() ); 
	  // std::cout << SyncOut() << std::this_thread::get_id() << " adding " << getColor()  << ","<< i << std::endl; 
	  for(nat j = 0; j < hisData.size(); ++j)
	    result[j] += hisData[j]; 
	}

      auto wasRead = msg.clearMessage(1); 
      if(not wasRead)
	{
	  std::cout << SyncOut() << "problem with color "<< getColor()  << "," << myRank << SHOW(size()) << std::endl; 
	  assert(wasRead ); 
	}

    }
  else 
    {
      while(not msg.clearMessage(1))
	; 
    }

  return result; 
}



template<typename T> 
auto LocalComm::gatherVariableLength(std::vector<T> myData, int root )  
  ->std::vector<T> 
{
  int myIdx = getIdx();
  int myCol = _colors.at(myIdx);
  auto result = std::vector<T>(); 
  auto& myMsg = _messages.at(myIdx); 

  myMsg.post(myData); 
  auto myGen = myMsg.getGeneration() ; 

  if(getRank() == root)
    {
      for(int i = 0; i < size() ;++i)
	{
	  int hisIdx = getIdx(myCol, i); 
	  auto &hisMsg = _messages.at(hisIdx); 
	  while(myGen != hisMsg.getGeneration())
	    ; 
	  auto msg = hisMsg.getMessage<T>();
	  result.insert(end(result), begin(msg), end(msg) ); 
	}

      auto wasCleared = myMsg.clearMessage(1); 
      assert(wasCleared); 

    }
  else 
    {
      while(not myMsg.clearMessage(1))
	; 
    }

  return result; 
}


template<typename T>
std::vector<T> LocalComm::allReduce(std::vector<T> myData)
{
  int root = 0 ; 
  myData = reduce(myData, root);
  myData = broadcast(myData,root); 
  return myData; 
}


template<typename T>
void LocalComm::postAsyncMessage(const std::vector<T> &message, int numRead, int tag)
{
  _asyncMessages.post(message, numRead, tag);
}


template<typename T>
std::tuple<bool,std::vector<T> > LocalComm::readAsyncMessage(int tag)
{
  return _asyncMessages.consumeOne<T>(tag);
}
