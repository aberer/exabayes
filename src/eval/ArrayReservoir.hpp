#ifndef REVERVE_ARRAYS_HPP
#define REVERVE_ARRAYS_HPP

#include <memory>
#include <map>
#include <unordered_map>
#include <vector>

#include "common.h"


typedef void* array_reservoir_t ; 
extern "C" 
{
  double* allocate(array_reservoir_t self, size_t a ); 
  void deallocate(array_reservoir_t self, double *array); 
}



class ArrayReservoir
{
public: 
  ArrayReservoir(bool freeOlds) : _freeOlds(freeOlds){}
  ~ArrayReservoir(); 

  ArrayReservoir(const ArrayReservoir& rhs) = delete; 
  ArrayReservoir& operator=(const ArrayReservoir& rhs) = delete; 

  double* allocate(size_t requiredLength ); 
  void deallocate(double* array); 

private: 
  std::unordered_map<double*,nat> _usedArrays; 
  std::multimap<nat,double*> _unusedArrays; 
  
  static const double numGammaCats; 
  static const double thresholdForNewSEVArray; 
  bool _freeOlds; 
}; 
#endif