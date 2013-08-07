#include <iostream>
#include <cassert>
#include "SuccessCounter.hpp"

#include "common.h"

SuccessCounter::SuccessCounter() 
  : globalAcc(0)
  , globalRej(0)
  , localAcc(0)
  , localRej(0)
  , batch(0)
{
}


SuccessCounter::SuccessCounter(const SuccessCounter& rhs)
  : globalAcc(rhs.globalAcc)
  , globalRej(rhs.globalRej)
  , localAcc(rhs.localAcc)
  , localRej(rhs.localRej)
  , batch(rhs.batch)
{  
}


SuccessCounter SuccessCounter::operator=( const SuccessCounter &rhs) 
{
  if(&rhs == this)
    return *this; 
  else 
    return SuccessCounter(rhs);   
} 



void SuccessCounter::accept() 
{
  globalAcc++;
  localAcc++; 
  addOrReplace(true); 
}


void SuccessCounter::reject()
{
  globalRej++;
  localRej++; 
  addOrReplace(false); 
}


double SuccessCounter::getRatioInLast100() const 
{
  double acc = 0, rej = 0 ; 
  for(auto b : lastX )
    if(b)
      acc++;
    else 
      rej++; 
  return acc / (acc+rej) ; 
}



double SuccessCounter::getRatioInLastInterval() const 
{ 
  double ratio =   ((double)localAcc) / ((double)(localAcc + localRej) ); 
  // cout << "ratio is " <<  ratio << endl; 
  return  ratio ; 
}


double SuccessCounter::getRatioOverall() const 
{
  return (double) globalAcc / ((double)(globalAcc + globalRej)); 
}

void SuccessCounter::nextBatch()
{
  batch++;
  reset();
}



void SuccessCounter::reset()
{
  localRej = 0; localAcc = 0; 
}


void SuccessCounter::addOrReplace(bool acc)
{ 
  if(lastX.size() ==  SIZE_OF_LAST) 
    {
      lastX.pop_front(); 
    }
  lastX.push_back(acc); 
}


void SuccessCounter::readFromCheckpoint( std::istream &in )   
{
  globalAcc = cRead<int>(in); 
  globalRej = cRead<int>(in); 
  localAcc = cRead<int>(in); 
  localRej = cRead<int>(in); 
  batch = cRead<int>(in); 
}
 
void SuccessCounter::writeToCheckpoint( std::ostream &out)  const
{
  cWrite<int>(out, globalAcc); 
  cWrite<int>(out, globalRej); 
  cWrite<int>(out, localAcc); 
  cWrite<int>(out, localRej); 
  cWrite<int>(out, batch); 
} 



SuccessCounter SuccessCounter::operator+(const SuccessCounter &rhs) const 
{
#ifdef UNSURE
  // this whole addition thing is not really consistent. It does its job for swap matrices, I guess. 
  assert(0); 
#endif

  SuccessCounter result; 
  result.globalAcc = rhs.globalAcc + globalAcc; 
  result.globalRej = rhs.globalRej + globalRej; 
  result.localAcc = rhs.localAcc + localAcc; 
  result.localRej = rhs.localRej + localRej; 
  result.batch = rhs.batch + batch; 

  return result; 
} 



