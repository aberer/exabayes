#ifndef _MY_EXTENSIONS_HPP
#define _MY_EXTENSIONS_HPP

#include <vector>
#include <type_traits>
#include <memory>
#include <iostream>

#include "log_double.hpp"

#include "common.h"

#define PRINT(x) (tout << "#x=x"  << std::endl) 


// still missing in c++11. 
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}


// a bit more number formatting  
template<typename T> class ThousandsSeparator : public std::numpunct<T> 
{
public:
  ThousandsSeparator(T Separator) : m_Separator(Separator) {}

protected:
  std::string do_grouping() const
  {
    return "\03";
  }

private:
    T m_Separator;
}; 

void formatRange(std::ostream &out, const std::vector<nat> &values) ; 

template<typename T> struct shared_pod_ptr
{
  shared_pod_ptr(T *arg = nullptr)
    : _impl{arg, [](T* bla){ free(bla) ; }}
  {}
  
  T* get() {return _impl.get() ; }

private: 
  std::shared_ptr<T> _impl; 
}; 

template<typename T, int ALIGN> T* aligned_malloc( size_t size ); 



#include "AlignedAllocator.hpp"


template<typename T>
struct aligned_vector
{
typedef std::vector<T, AlignedAllocator<T, EXA_ALIGN> >  type; 
}; 

// that's mostly for development. accessing an environment variable is easier than 
std::string getEnvironmentVariable( const std::string  & key ) ; 




// okay, so doubles are not allowed as a non-type template parameter =( 

// template<typename T, T &MIN_VAL, T &MAX_VAL, bool MIN_OPEN = false, bool MAX_OPEN = false >
// class is_in_range
// {
// public: 
//   is_in_range(T num) 
//     : _num{num}
//   {
//   }

//   operator bool()
//   {
//     auto result = true; 
    
//     result &= MIN_OPEN ? (MIN_VAL < _num) : (MIN_VAL <= _num ) ; 
//     result &= MAX_OPEN ? (_num < MAX_VAL) : (_num <= MAX_VAL); 

//     return result; 
//   }


// private: 
//   T _num; 
  
// };


#endif
