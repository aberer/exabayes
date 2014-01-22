#ifndef _MY_EXTENSIONS_HPP
#define _MY_EXTENSIONS_HPP

#include <memory>

// still missing in c++11. 

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

#endif
