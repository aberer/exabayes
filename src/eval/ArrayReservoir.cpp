#include "eval/ArrayReservoir.hpp" 
#include "GlobalVariables.hpp" 
#include <algorithm>
#include <cassert>
#include "axml.h"

const double ArrayReservoir::numGammaCats = 4; 
const double ArrayReservoir::thresholdForNewSEVArray = 1.05 ; 


// adding c-bindings
double* allocate(array_reservoir_t self, size_t required)
{
  auto castedSelf = static_cast<ArrayReservoir*>(self); 
  return castedSelf->allocate(required); 
}


void deallocate(array_reservoir_t self, double *array)
{
  auto castedSelf = static_cast<ArrayReservoir*>(self); 
  castedSelf->deallocate(array);
}


ArrayReservoir::~ArrayReservoir()
{
  for(auto elem : _unusedArrays)
    exa_free(std::get<1>(elem));
} 



// nat gctr = 0; 

double* ArrayReservoir::allocate(size_t requiredLength)
{
  double *result = nullptr; 

  auto found = _unusedArrays.lower_bound(requiredLength); 
  
  auto lengthOfFound = std::get<0>(*found);
  
  // we free anyway such that we do not acccumulate a lot of large
  // arrays when the SEV technique is used
  bool freeOldArrayDespiteAvailable =  _freeOlds &&   ( requiredLength * thresholdForNewSEVArray  < lengthOfFound )  ; 

  if( found != end(_unusedArrays)
      && not freeOldArrayDespiteAvailable )
    {
      auto elem = *found; 
      _usedArrays.insert(std::make_pair(std::get<1>(elem), std::get<0>(elem)));
      _unusedArrays.erase(found);
      result = std::get<1>(elem); 
      // tout  << _unusedArrays.size() + _usedArrays.size() << "=" << _unusedArrays.size() << "/" << _usedArrays.size()  << " => REUSE "  << std::endl; 
    }
  else 
    {				
      // tout  << _unusedArrays.size() + _usedArrays.size() << "=" << _unusedArrays.size() << "/" << _usedArrays.size()  << " => NEW "  << std::endl; 
      
      // if(freeOldArrayDespiteAvailable)
      // 	{
	  // ++gctr; 
	  // tout <<  gctr << "\tfreeing an array anyway. Needed "  << requiredLength << ", but smalles was " << lengthOfFound << std::endl; 
	// }

      if(_freeOlds &&   _unusedArrays.size() > 0  )
	{
	  auto iter =  _unusedArrays.end();
	  --iter; 
	  exa_free(std::get<1>(*iter));
	  _unusedArrays.erase(iter); 
	}

      result = (double*) exa_malloc_aligned(requiredLength);
      _usedArrays.insert(std::make_pair(result, requiredLength));
    }

  assert(result != nullptr) ; 
  return result; 
}


void ArrayReservoir::deallocate(double* array)
{
  assert(array != nullptr); 

  auto iter = _usedArrays.find(array); 

  assert(iter != _usedArrays.end()); 
  auto elem = *iter; 

  _unusedArrays.insert(std::make_pair(std::get<1>(elem), std::get<0>(elem))); 
  _usedArrays.erase(iter);

  assert(_usedArrays.find(array) == _usedArrays.end() ); 

  // tout << _unusedArrays.size() + _usedArrays.size() << "=" << _unusedArrays.size() << "/" << _usedArrays.size()  << " => PUT BACK "  << std::endl; 

  array = nullptr; 
} 
