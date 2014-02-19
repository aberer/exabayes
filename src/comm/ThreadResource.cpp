#include "comm/ThreadResource.hpp"

#include "comm/ParallelSetup.hpp"
#include "config/CommandLine.hpp"
#include <cstring>

#include <sched.h>

void exa_main ( CommandLine &cl,  ParallelSetup* pl ); 
 


ThreadResource::ThreadResource(CommandLine& cl, ParallelSetup* pl)
  : _threadsAreReleased(false)
{
  // this starts the threads and makes them wait in the function threadStart
  auto foo = std::bind(&ThreadResource::threadStart,  this, std::ref(cl), pl);
  int highestRank = 0; 
  _tid2rank[MY_TID] = highestRank++; 
  for(int i = 1 ; i < cl.getNumThreads()  ; ++i)
    {
      _threads.push_back(std::thread(foo  )); 
      _tid2rank[_threads.back().get_id()] = highestRank++; 
    }
  
}



void swap(ThreadResource &lhs, ThreadResource& rhs)
{
  using std::swap; 
  swap(lhs._threads, rhs._threads); 
  swap(lhs._tid2rank, rhs._tid2rank); 
  swap(lhs._threadsAreReleased, rhs._threadsAreReleased); 
} 


ThreadResource& ThreadResource::operator=(ThreadResource rhs) 
{
  swap(*this, rhs); 
  return *this; 
} 


void ThreadResource::threadStart(CommandLine& cl, ParallelSetup* pl)
{
  // auto off = pl->getGlobalComm().getOffsetForThreadPin() * _threads.size(); 
  // setAffinity(off);

  while(not _threadsAreReleased)
    ; 

  exa_main(cl, pl); 
}


void ThreadResource::releaseThreads()
{
  _threadsAreReleased = true; 
}


ThreadResource::~ThreadResource()
{
  nat ctr = 0; 
  for(auto &t : _threads)
    t.join();
}


void ThreadResource::setAffinity(int offset)
{
  cpu_set_t  mask;
  memset(&mask, 0, sizeof(cpu_set_t)); 

  CPU_SET(_tid2rank[MY_TID] + offset, &mask);
  int result = sched_setaffinity(0, sizeof(mask), &mask); 

  assert( ( result == 0)   && "could not set process affinity");
}
