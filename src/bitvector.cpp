#include "bitvector.hpp"

#include <iostream>
#include <cassert>

using std::endl;
using std::cout; 

const bv_elem bitvector::static_mask[32] =
  {
    1u << 0,      1u << 1,      1u << 2,      1u << 3,
    1u << 4,      1u << 5,      1u << 6,      1u << 7,
    1u << 8,      1u << 9,      1u << 10,      1u << 11,
    1u << 12,      1u << 13,      1u << 14,      1u << 15,
    1u << 16,      1u << 17,      1u << 18,      1u << 19,
    1u << 20,      1u << 21,      1u << 22,      1u << 23,
    1u << 24,      1u << 25,      1u << 26,      1u << 27,
    1u << 28,      1u << 29,      1u << 30,      1u << 31
  };


const bv_elem bitvector::allOne =   ( ( (1 << 16) -  1 )  << 16  ) +  ((1 << 16) - 1 ) ; 
const bv_elem bitvector::allZero = 0u; 


const bv_elem bitvector::nOne[32] =
  {
    0u, 
    bitvector::allOne << 31, bitvector::allOne << 30, bitvector::allOne << 29, bitvector::allOne << 28,
    bitvector::allOne << 27, bitvector::allOne << 26, bitvector::allOne << 25, bitvector::allOne << 24,
    bitvector::allOne << 23, bitvector::allOne << 22, bitvector::allOne << 21, bitvector::allOne << 20,
    bitvector::allOne << 19, bitvector::allOne << 18, bitvector::allOne << 17, bitvector::allOne << 16,
    bitvector::allOne << 15, bitvector::allOne << 14, bitvector::allOne << 13, bitvector::allOne << 12,
    bitvector::allOne << 11, bitvector::allOne << 10, bitvector::allOne << 9, bitvector::allOne << 8,
    bitvector::allOne << 7, bitvector::allOne << 6, bitvector::allOne << 5, bitvector::allOne << 4,
    bitvector::allOne << 3, bitvector::allOne << 2, bitvector::allOne << 1
  };



void bitvector::resize(size_t theSize)
{
  auto needed = elemSize(theSize);
  // std::cout << "resizing to " << needed  << std::endl;
  _array.resize(needed);
  _usedSize = theSize;
}


size_t bitvector::elemSize(size_t theSize) 
{
  return theSize / bitsPerElem +   ( ( theSize % bitsPerElem  ) == 0 ? 0 : 1  ) ;
}


size_t bitvector::size() const
{
  return _usedSize; 
}


bool bitvector::get(size_t index) const
{
  if( size() <= index)
    return false; 
  return _array[index / bitsPerElem ] & static_mask[ index % bitsPerElem]; 
}


void bitvector::ensureCapacity(size_t c)
{
  auto index = c ; 
  if(_array.size() * bitsPerElem <= index )
    resize( index );
  else if( _usedSize < c  )
    _usedSize = c; 
}

void bitvector::set(size_t index ) 
{
  ensureCapacity(index + 1 );
  _array.at( index / bitsPerElem) |= static_mask[ index % bitsPerElem ];
}

void bitvector::unset(size_t index)
{
  ensureCapacity(index + 1 );
  _array[index / bitsPerElem] &=  (~ static_mask[ index % bitsPerElem ]); 
}


size_t bitvector::count() const
{
  // TODO 
  auto result = 0u;

  for(auto i = 0u; i < _usedSize; ++i )
    {
      if(get(i))
        ++result;
    }
  
  return result; 
}



bitvector operator|( bitvector const& lhs, bitvector const& rhs)
{
  if(&rhs == &lhs)
    {
      cout << "\n\ndirect return\n\n" << std::endl; 
      return lhs;
    }

  auto *larger = &lhs;  
  auto *other =  &rhs; 
  if( larger->size() < other->size() )
    std::swap(larger, other); 

  auto result =  *larger; 

  auto resultIter = result._array.begin();
  auto resultOther = other->_array.begin();

  while(resultOther != other->_array.end())
    {
      *resultIter |= *resultOther;
      ++resultIter;
      ++resultOther; 
    }

  return result; 
}


bitvector operator&(  bitvector const& lhs, bitvector const& rhs)
{
  auto result = bitvector();
  auto maxSize = std::max( lhs.size() , rhs.size());
  result.resize( maxSize ); 

  
  auto minSize = std::min(lhs.size(), rhs.size()) ;

  auto minElemSize = bitvector::elemSize(minSize); 

  assert(lhs._array.size() >= minElemSize ) ;
  assert(rhs._array.size() >= minElemSize ) ; 
  
  for(auto i =0u ; i < minElemSize; ++i)
    result._array.at(i) = lhs._array.at(i) & rhs._array.at(i);

  return result;
}


bitvector operator-(bitvector const& lhs, bitvector const rhs)
{
  return lhs & (~ rhs) ; 
}


bitvector& bitvector::operator|=(bitvector const& rhs )
{
  *this = *this | rhs; 
  return *this; 
}

bitvector& bitvector::operator&=(bitvector const& rhs )
{
  *this = *this & rhs ; 
  return *this; 
}


std::ostream& operator<<(std::ostream& s, bitvector const& rhs )
{
  for(auto i = 0u; i < rhs.size() ; ++i)
    {
    if(rhs.get(i))
      s << "1" ;
    else
      s << "0"; 
    }
  
  return s; 
}


void bitvector::iterator::set(bool val)
{
  if(val)
    _ref->set(_index) ;
  else
    _ref->unset(_index) ;
}


bool bitvector::iterator::get() const
{
  return _ref->get(_index);
}


bitvector bitvector::operator~() const 
{
  auto cpy = *this; 
  for(auto &v : cpy._array)
    v = ~v; 
  return cpy; 
}


bool bitvector::operator==( bitvector const&other) const
{
  auto minSize = elemSize(std::min(size(), other.size()));
  auto maxSize = elemSize(std::max(size(), other.size()));

  for(auto i = 0u; i < minSize; ++i)
    {
      if( _array[i] != other._array[i] )
        return false; 
    }

  auto &larger =  size() > other.size() ? *this : other;
  for(auto i = minSize; i < maxSize; ++i)
    {
      if(larger._array[i] != 0 )
        return false;
    }

  return true; 
}


bitvector::bitvector( )
  :_array{}
  , _usedSize{0u}
{
}





bitvector::bitvector( std::initializer_list<int> l )
  : bitvector()
{
  for(auto v : l )
    set(v); 
}




bool bitvector::operator<(bitvector const& rhs) const
{
  return  ( ( *this & rhs) == *this )  && count() < rhs.count() ; 
}


bitvector::operator vector<size_t>() const 
{
  auto result = vector<size_t>{};

  for(auto i = 0u ;i < size() ; ++i)
    {
      if( get(i))
        result.push_back(i);
    }
  
  return result; 
}


bitvector bitvector::symmetricDifference( bitvector const& rhs) const
{
  return (*this | rhs ) - (*this & rhs); 
}
