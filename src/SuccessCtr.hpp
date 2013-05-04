#ifndef _SUCCESSCTR_H
#define _SUCCESSCTR_H

#include <fstream>
#include <iomanip>
#include <list>

using namespace std; 


#define SIZE_OF_LAST 100 

// TODO  an success interval would be nice 


class SuccessCtr
{
public: 
  SuccessCtr() : globalAcc(0),globalRej(0), localAcc(0), localRej(0), batch(0){}

  void accept(){globalAcc++;localAcc++; addOrReplace(true); }
  void reject(){globalRej++;localRej++; addOrReplace(false); }

  int getRecentlySeen(){return localAcc + localRej; }
  

  double getRatioInLast100() const 
  {
    double acc = 0, rej = 0 ; 
    for(auto b : lastX )
      if(b)
	acc++;
      else 
	rej++; 
    return acc / (acc+rej) ; 
  }
  double getRatioInLastInterval() const  { return (double)localAcc / ((double)(localAcc + localRej) ); }
  double getRatioOverall() const {return (double) globalAcc / ((double)(globalAcc + globalRej)); }

  void nextBatch(){batch++;reset();}
  int getBatch(){return batch; }

  friend ostream& operator<<(ostream& rhs, const SuccessCtr &b ); 

private: 
  void reset(){localRej = 0; localAcc = 0; }
  void addOrReplace(bool acc)
  { 
    if(lastX.size() ==  SIZE_OF_LAST) 
      {
	lastX.pop_front(); 
      }
    lastX.push_back(acc); 
  }
  
  
  
  list<bool> lastX ; 		// last x events. Only for debug, no functionality in tuning 

  int globalAcc; 
  int globalRej; 
  
  int localAcc; 
  int localRej; 
  int batch; 
}; 



inline ostream& operator<<(ostream& rhs, const SuccessCtr &sctr)
{  
  return  rhs <<  setprecision(1) <<  fixed << sctr.globalAcc  << "/" << sctr.globalRej << " (" <<  sctr.getRatioInLast100() * 100 << "%/"  << sctr.getRatioOverall() *  100    << "%)" ;  
}

#endif
