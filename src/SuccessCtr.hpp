#ifndef _SUCCESSCTR_H
#define _SUCCESSCTR_H

#include <fstream>
#include <iomanip>

using namespace std; 


// TODO  an success interval would be nice 


class SuccessCtr
{
public: 
  SuccessCtr() : globalAcc(0),globalRej(0), localAcc(0), localRej(0), batch(1){}

  void accept(){globalAcc++;localAcc++;}
  void reject(){globalRej++;localRej++;}

  int getRecentlySeen(){return localAcc + localRej; }
  
  double getRatioInLastInterval() const  { return (double)localAcc / ((double)(localAcc + localRej) ); }
  double getRatioOverall() const {return (double) globalAcc / ((double)(globalAcc + globalRej)); }

  void nextBatch(){batch++;reset();}
  int getBatch(){return batch; }

  void reset(){localRej = 0; localAcc = 0; }
  
  friend ostream& operator<<(ostream& rhs, const SuccessCtr &b ); 

private: 
  int globalAcc; 
  int globalRej; 
  
  int localAcc; 
  int localRej; 
  int batch; 
}; 



inline ostream& operator<<(ostream& rhs, const SuccessCtr &sctr)
{  
  return rhs << sctr.globalAcc  << "/" << sctr.globalRej << " (" <<  setprecision(1) << fixed << sctr.getRatioOverall() *  100    << "%)" ;  
}

#endif
