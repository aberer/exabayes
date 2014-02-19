#ifndef _LOCAL_COMM_HPP
#define _LOCAL_COMM_HPP

#include <unordered_map>
#include "comm/Message.hpp"
#include "comm/threadDefs.hpp"
#include <iostream>

#include "comm/MessageQueue.hpp"

class LocalComm
{
  typedef LocalComm SELF; 

public: 
  /** 
      @param tid2rank contains absolute ranks of all threads 
   */ 
  LocalComm(std::unordered_map<tid_t,int> tid2rank); 
  LocalComm(const LocalComm& rhs); 
  LocalComm(LocalComm &&rhs) ; 
  LocalComm& operator=(LocalComm rhs ); 

  friend void swap(LocalComm& lhs, LocalComm& rhs ); 

  friend std::ostream& operator<<(std::ostream& out, const LocalComm& rhs); 

  std::unordered_map<tid_t,int> getTid2Ranking  () const {return _tid2LocCommIdx; }
  void setColors(std::vector<int> colors) { _colors = colors; }
  void setRanks(std::vector<int> ranks) {_ranks = ranks; }
  
  #include "comm/CommCore.hpp"

  int getIdx() const {return _tid2LocCommIdx.at(MY_TID); }
  int getNumThreads() const ; 

  int getIdx(int col, int rank) const ; 
  int getColor() const {return _colors.at(_tid2LocCommIdx.at(MY_TID)); }

  bool checkAsyncMessage( int tag )  const ; 

  template<typename T>
  void postAsyncMessage(const std::vector<T> &message, int numRead, int tag); 
  template<typename T>
  std::tuple<bool,std::vector<T> > readAsyncMessage(int tag); 

private: 			// ATTRIBUTES
  std::vector<Message> _messages; 
  std::unordered_map<tid_t,int> _tid2LocCommIdx;

  std::vector<int> _colors; 
  std::vector<int> _ranks; 
  int _size; 

  MessageQueue _asyncMessages; 
}; 


#include "comm/LocalComm.tpp"

#endif
