#ifndef _MY_EXTENSIONS_HPP
#define _MY_EXTENSIONS_HPP

#include <vector>
#include <memory>
#include <iomanip>

#include <iostream>

#include "common.h"

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
  // T do_thousands_sep() const  {
  //   return m_Separator;
  // }
  
  std::string do_grouping() const
  {
    return "\03";
  }

private:
    T m_Separator;
}; 



void formatRange(std::ostream &out, const std::vector<nat> &values) ; 
// {
//   bool inRange = false; 
//   for(auto iter = begin(values); iter != end(values) ; ++iter )
//     {
//       if(iter == begin(values))
// 	out << *iter; 
//       else if(iter == end(values) -1 )
// 	{
// 	  if(inRange)
// 	    out << "-" << *iter ;  
// 	  else 
// 	    out << "," << *iter; 
// 	}
//       else 
// 	{ 
// 	  bool haveBeenInRange = ((iter-1) ==  begin(values)) || inRange; 
// 	  inRange = *iter == *(iter-1) + 1;  

// 	  if(not haveBeenInRange )
// 	    {
// 	      out << "," << *iter ; 
// 	    }
// 	  else
// 	    {
// 	      if(not inRange)
// 		out << "-" << *iter; 
// 	    }
// 	}
//     }
  
// }


#endif
