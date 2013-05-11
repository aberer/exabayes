#ifndef _SUCCESSCTR_H
#define _SUCCESSCTR_H

#include <fstream>
#include <iomanip>
#include <list>

using namespace std; 


#define SIZE_OF_LAST 100 

// TODO  an success interval would be nice 


class SuccessCounter
{
public: 
  SuccessCounter();
  SuccessCounter(const SuccessCounter& rhs); 

  void accept(); 
  void reject();
  int getRecentlySeen(){return localAcc + localRej; }
  double getRatioInLast100() const ; 
  double getRatioInLastInterval() const ; 
  double getRatioOverall() const ; 
  void nextBatch(); 
  int getBatch(){return batch; }

  friend ostream& operator<<(ostream& rhs, const SuccessCounter &b ); 

private:   
  list<bool> lastX ; 		// last x events. Only for debug, no functionality in tuning 

  int globalAcc; 
  int globalRej; 
  
  int localAcc; 
  int localRej; 
  int batch; 

  void reset();  
  void addOrReplace(bool acc); 
}; 


inline ostream& operator<<(ostream& rhs, const SuccessCounter &sctr)
{  
  return  rhs <<  setprecision(1) <<  fixed << sctr.globalAcc  << "/" << sctr.globalRej << " (" <<  sctr.getRatioInLast100() * 100 << "%/"  << sctr.getRatioOverall() *  100    << "%)" ;  
}

#endif
