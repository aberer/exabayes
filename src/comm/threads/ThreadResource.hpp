#ifndef _THREAD_RESOURCE_HPP
#define _THREAD_RESOURCE_HPP

#include "comm/threads/threadDefs.hpp"

#include <memory>
#include <vector>
#include <unordered_map>

#include "common.h"

class CommandLine; 
class ParallelSetup; 

class ThreadResource
{
public: 
  ThreadResource(CommandLine &cl, ParallelSetup* pl); 
  ThreadResource(ThreadResource&& rhs) = default; 
  ThreadResource& operator=(ThreadResource rhs) ; 

  friend void swap(ThreadResource &lhs, ThreadResource& rhs); 
  ~ThreadResource(); 

  int getNumThreads() const {return _threads.size() + 1; }

  void threadStart(CommandLine& cl, ParallelSetup* pl); 
  void releaseThreads(); 
  std::unordered_map<std::thread::id,int>getTid2Ranking() const {return _tid2rank; } 
  
  void pinThreads(int numProcPerNode, int remoteRank); 

private: 
  std::vector<std::thread> _threads; 
  std::unordered_map<tid_t, int> _tid2rank; // thread rank 
  volatile bool _threadsAreReleased; 
}; 

#endif
