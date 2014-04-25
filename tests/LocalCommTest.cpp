#include "comm/LocalComm.hpp"

#include <vector>
#include <thread>
#include <unordered_map>


struct fence
{
  LocalComm *comm; 
  std::mutex mtx; 
  bool isThere; 
  std::function<void(LocalComm&)> fun; 
}; 



void start(fence &fnc)
{
  std::lock_guard<std::mutex>(fnc.mtx);
		  fnc.isThere = true; 
}

bool isThere(fence &fnc)
{
  std::lock_guard<std::mutex>(fnc.mtx);
		  return fnc.isThere; 
}



void execute(  fence *myFencePtr)
{
  fence &myFence = *myFencePtr; 
  while(not isThere(myFence)  )
    ; 

  auto &comm = *myFence.comm; 
  myFencePtr->fun(comm); 
}




void executeWithThreads( int numThreads, std::function<void(LocalComm&)> lam)
{
  auto threads = std::vector<std::thread>{}; 
  
  fence* myFence = new fence;  
  myFence->comm = nullptr; 
  myFence->isThere = false; 
  auto&& fun = std::bind(execute, myFence); 

  auto ids = std::unordered_map<std::thread::id,int>{}; 

  ids.insert(std::make_pair(MY_TID, 0)); 
  for(int i = 1; i < numThreads ; ++i)
    {
      threads.push_back(std::thread(fun)); 
      ids.insert(std::make_pair(threads.back().get_id(), i)); 
    }

  auto&& comm = LocalComm(ids); 
  myFence->comm = &comm; 
  myFence->fun = lam; 
  assert(comm.isValid()); 
  start(*myFence);
  
  execute(myFence);

  for(auto &t : threads)
    t.join();
}



void executeWithThreadsAfterSplit( int numThreads, std::function<void(LocalComm&)> lam, std::vector<int> color, std::vector<int> rank)
{
  auto threads = std::vector<std::thread>{}; 
  
  fence* myFence = new fence;  
  myFence->comm = nullptr; 
  myFence->isThere = false; 
  auto&& fun = std::bind(execute, myFence); 

  auto ids = std::unordered_map<std::thread::id,int>{}; 

  ids.insert(std::make_pair(MY_TID, 0)); 
  for(int i = 1; i < numThreads; ++i)
    {
      threads.push_back(std::thread(fun)); 
      ids.insert(std::make_pair(threads.back().get_id(), i)); 
    }

  auto&& comm = LocalComm(ids); 
  comm = comm.split(color, rank);
  myFence->comm = &comm; 
  myFence->fun = lam; 
  assert(comm.isValid()); 
  start(*myFence);
  
  execute(myFence);

  for(auto &t : threads)
    t.join();
}



TEST(LocalCommTest, bcastTest)
{
  auto lam = [](LocalComm& comm) 
    {
      auto myData = std::vector<nat>(10,0); 
      if(comm.getRank() == 0)
	{
	  for(nat i =0; i < myData.size();  ++i)
	    myData[i] = i+1; 
	}

      myData = comm.broadcast(myData,0); 

      for(nat i = 0; i< myData.size() ;++i)
	ASSERT_EQ(i+1, myData[i] ); 
    }; 

  for(int i = 2; i < 9 ; ++i)
    {
      executeWithThreads(i , lam); 
      // std::cout << SyncOut() << "================ okay ================ "<< std::endl; 
    }
}




TEST(LocalCommTest, bcastTestSplitted)
{
  auto lam = [](LocalComm& comm) 
    {
      auto myData = std::vector<nat>(10,0); 
      if(comm.getRank() == 0 )
	{
	  for(nat i =0 ; i< myData.size(); ++i)
	    myData[i] = i + comm.getColor(); 
	}

      myData = comm.broadcast(myData, 0);

      for(nat i = 0; i < myData.size(); ++i)
	ASSERT_TRUE(myData[i] == i + comm.getColor()); 
    }; 
  
  auto cols = std::vector<int>{0,0,1,1}; 
  auto ranks = std::vector<int>{0,1,0,1}; 

  executeWithThreadsAfterSplit(4, lam, cols,ranks);
}



TEST(LocalCommTest, scatterv)
{
  auto lam = [](LocalComm& comm)  
    {
      int ROOT = 0; 
      auto counts = std::vector<int>{3,2,3,2};
      auto displ = std::vector<int>{0,3,5,8}; 
      auto myData = std::vector<int>(counts.at(comm.getRank()),0); 
      auto theData = std::vector<int>();       
      if(comm.getRank() == ROOT) 
	{

	  nat ctr = 0; 
	  for(auto &c : counts)
	    {
	      for(int i = 0; i < c ; ++i)
		theData.push_back(ctr); 
	      ++ctr; 
	    }
	  // std::cout << "the data is " << theData << std::endl; 
	  assert(theData.size() == 10); 
	}
      assert(comm.size() == 4); 

      myData = comm.scatterVariableKnownLength(theData,  counts,displ, ROOT );
      
      for(auto &d : myData )
	{
	  int myRank = comm.getRank();  
	  ASSERT_EQ(d, myRank);
	}
    }; 
  
  executeWithThreads(4, lam);
}


TEST(LocalCommTest, splitScatterv)
{
  auto lam = [](LocalComm& comm)  
    {
      if(comm.getColor() == 0 )
	return; 

      int ROOT = 0; 
      auto counts = std::vector<int>{3,2,3,2};
      auto displ = std::vector<int>{0,3,5,8}; 
      auto myData = std::vector<int>(counts.at(comm.getRank()),0); 
      auto theData = std::vector<int>();       
      if(comm.getRank() == ROOT) 
	{
	  nat ctr = 0; 
	  for(auto &c : counts)
	    {
	      for(int i = 0; i < c ; ++i)
		theData.push_back(ctr); 
	      ++ctr; 
	    }
	  // std::cout << "the data is " << theData << std::endl; 
	  assert(theData.size() == 10); 
	}
      assert(comm.size() == 4); 

      myData = comm.scatterVariableKnownLength(theData,  counts,displ, ROOT );
      
      for(auto &d : myData )
	{
	  int myRank = comm.getRank();  
	  ASSERT_EQ(d, myRank);
	}
    }; 
  
  auto cols = std::vector<int>{0,0,0,0,1,1,1,1}; 
  auto ranks = std::vector<int>{0,1,2,3,0,1,2,3};

  executeWithThreadsAfterSplit(8, lam, cols, ranks); 
}

static int rankSum(int size)
{
  auto result = 0; 
  for(int i = 0; i< size; ++i)
    result += i; 
  return result; 
}

TEST(LocalCommTest, reduceTest)
{
  auto lam = [](LocalComm& comm)
    {
      int ROOT = 0; 
      auto myData = std::vector<int>{comm.getRank(), comm.getRank() + 1 };
      myData = comm.reduce(myData,ROOT); 
      
      if(comm.getRank() == ROOT)
	{
	  ASSERT_EQ(rankSum(comm.size()), myData[0] ); 
	  ASSERT_EQ(rankSum(comm.size()) + comm.size(), myData[1] ); 
	}
    }; 
  
  executeWithThreads(4, lam);
}



TEST(LocalCommTest, allReduceTest)
{
  auto lm = [](LocalComm& comm)
    {
      for(int i = 0; i < 5; ++i)
	{
	  auto myData = std::vector<int>{comm.getRank(), comm.getRank() + 1 };
      
	  myData = comm.allReduce(myData); 

	  ASSERT_EQ(rankSum(comm.size()), myData[0] ); 
	  ASSERT_EQ(rankSum(comm.size()) + comm.size(), myData[1] ); 
	}
    }; 

  for(int i = 0; i < 12; ++i)
    {
      executeWithThreads(i, lm);
    }
}



TEST(LocalCommTest, reduceTestSplit)
{
  auto lam = [](LocalComm& comm)
    {
      int ROOT = 0; 
      auto myData = std::vector<int>{comm.getRank(), comm.getRank() + 1 };
      auto origData = myData;  
      
      myData = comm.reduce(myData,ROOT); 
      
      if(comm.getRank() == ROOT)
      	{
      	  ASSERT_EQ(rankSum(comm.size()), myData[0]); 
      	  ASSERT_EQ(rankSum(comm.size()) + comm.size(), myData[1] ); 
      	}
      
      
      // std::cout << SyncOut() << std::this_thread::get_id() << "," << comm.getColor() << "," << comm.getRank()  << " done" << std::endl; 

      myData = origData; 

      // ROOT = 1; 
      myData = comm.reduce(myData,ROOT); 
      
      if(comm.getRank() == ROOT)
      	{
      	  ASSERT_EQ(rankSum(comm.size()), myData[0] ); 
      	  ASSERT_EQ(rankSum(comm.size()) + comm.size(), myData[1] ); 
      	}
    }; 

  auto cols = std::vector<int>{0,0,0,0,1,1,1,1}; 
  auto ranks = std::vector<int>{0,1,2,3,0,1,2,3};  

  executeWithThreadsAfterSplit(8, lam, cols,ranks);
}



TEST(LocalCommTest, gathervTest)
{
  auto lam = [](LocalComm& comm)
    {
      if(comm.getColor() == 0)
	return;

      int ROOT = 0; 
      auto myData = std::vector<int>{comm.getRank()};
      if(comm.getRank() % 2 == 0)
	myData.push_back(comm.getRank()); 

      auto allData = comm.gatherVariableLength(myData, ROOT);
      if(comm.getRank() == ROOT)
	{
	  nat ctr = 0; 
	  for(int i = 0; i < comm.size() ;++i)
	    {
	      ASSERT_EQ(allData[ctr++] , i); 
	      if(i % 2 == 0)
		ASSERT_EQ(allData[ctr++], i); 
	    }
	}
    }; 

  auto cols = std::vector<int>{0,0,0,0,1,1,1,1}; 
  auto ranks = std::vector<int>{0,1,2,3,0,1,2,3};

  executeWithThreadsAfterSplit(8, lam, cols, ranks); 
}


TEST(LocalCommTest, gatherVariableKnownLength_Split_Test)
{

  auto lam = [](LocalComm&  comm)
    {
      int ROOT = 0; 
      auto myData = std::vector<int>{comm.getRank()}; 
      if(comm.getRank() % 2 == 0 )
	myData.push_back(comm.getRank() + 1 );
      
      auto countsPerProc = std::vector<int>();
      for(int i = 0; i < comm.size() ; ++i)
	{
	  if(i % 2 == 0 )
	    countsPerProc.push_back(2); 
	  else 
	    countsPerProc.push_back(1); 
	}

      auto displPerProc = std::vector<int>{0};
      for(int i =1; i < comm.size() ; ++i)
	displPerProc.push_back(displPerProc[i-1] + countsPerProc[i-1]); 

      auto allData = comm.gatherVariableKnownLength(myData, countsPerProc, displPerProc, ROOT);

      // if(comm.getColor() == 1 )
      // 	std::cout << SyncOut() << "i survived it! \n"; 

      if(comm.getRank() == ROOT)
	{
	  auto iter = begin(allData); 
	  nat ctr = 0; 
	  while(iter != end(allData))
	    {
	      ASSERT_EQ(ctr, *iter); 
	      ++iter; 
	      if(ctr % 2 == 0 )
		{
		  ASSERT_EQ(ctr+1, *iter); 
		  ++iter; 
		}
	      ++ctr;
	    }
	}

      // std::cout << SyncOut()  << comm.getColor() << "," << comm.getRank() << " finished" << std::endl; 

    }; 

  auto cols = std::vector<int>{0,0,0,0,1,1,1,1}; 
  auto ranks = std::vector<int>{0,1,2,3,0,1,2,3}; 
  executeWithThreadsAfterSplit(8,lam,cols,ranks); 
}
