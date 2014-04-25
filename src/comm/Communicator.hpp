#ifndef _NEW_COMMUNICATOR_HPP
#define _NEW_COMMUNICATOR_HPP

#include "comm/RemoteComm.hpp"      
#include "comm/LocalComm.hpp"
#include "comm/threads/threadDefs.hpp" 

#include <vector>
#include <cassert>

#include <unordered_map>

class ThreadResource; 

class Communicator
{
  typedef Communicator SELF ; 
public: 
  /** 
      @param tid2rank absolute ranks in the communicator ; contains only ranks of local threads   
  */ 
  Communicator(std::unordered_map<tid_t,int> tid2rank); 
  Communicator(const Communicator& rhs) = delete; 
  Communicator(Communicator &&rhs) = default; 
  ~Communicator(){}
  Communicator& operator=(Communicator rhs) ; 
  friend void swap(Communicator& lhs, Communicator &rhs); 

  void createSendRequest(std::vector<char> array, int dest, int tag, CommRequest& req);
  void createRecvRequest(int src, int tag, nat length, CommRequest& req); 

#include "comm/CommCore.hpp"

  friend std::ostream& operator<<(std::ostream & out, const Communicator& rhs); 
  
  int mapToLocalRank( int rank) const  ; 
  int mapToRemoteRank(int rank) const ; 

  int getProcsPerNode() ; 
  
  LocalComm&  getLocalComm() ; 
  RemoteComm& getRemoteComm(); 

  void initWithMaxChains(int numChains, int numThreadsChecking); 
  
private:
  RemoteComm _remoteComm; 
  LocalComm _localComm; 
};  


#include "comm/Communicator.tpp"

#endif


