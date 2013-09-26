#ifndef _SUCCESSCTR_H
#define _SUCCESSCTR_H

#include <iostream>
#include <iomanip>
#include <list>

typedef  unsigned int nat ; 

#include "Serializable.hpp"

#define SIZE_OF_LAST 100 

// TODO  an success interval would be nice 


class SuccessCounter : public Serializable
{
public: 
  SuccessCounter();  
  SuccessCounter(const SuccessCounter& rhs); 
  SuccessCounter operator=( const SuccessCounter &rhs) ; 

  void accept(); 
  void reject();
  int getRecentlySeen() const {return localAcc + localRej; }
  double getRatioInLastInterval() const ; 
  double getRatioOverall() const ; 
  void nextBatch(); 
  int getBatch() const {return batch; }
  nat getTotalSeen()const  {return globalAcc + globalRej; }
  
  virtual void deserialize( std::istream &in )   ; 
  virtual void serialize( std::ostream &out) const ;   


  SuccessCounter operator+(const SuccessCounter &rhs) const ; 
  
private: 			// METHODS
  void reset();  
  // void addOrReplace(bool acc); 

private:		  // ATTRIBUTES
  int globalAcc; 
  int globalRej;   
  int localAcc; 
  int localRej; 
  int batch; 			

  friend std::ostream& operator<<(std::ostream& rhs, const SuccessCounter &b ); 
}; 

#endif
