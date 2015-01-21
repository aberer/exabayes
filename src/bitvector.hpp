#ifndef BITVECTOR_H
#define BITVECTOR_H

#include <vector>
#include <cstddef>
#include <climits>
#include <exception>
#include <iosfwd>

using std::reference_wrapper; 
using std::vector;

using bv_elem = unsigned int ; 

class bitvector
{
public:
  class iterator; 
  
public: 
  static const bv_elem static_mask[32];
  static const bv_elem bitsPerElem = sizeof(int) *  CHAR_BIT;

  static const bv_elem allOne;
  static const bv_elem allZero;

  static const bv_elem nZero[32];
  static const bv_elem nOne[32];
  
public:
  bitvector( ); 
  bitvector( std::initializer_list<int> l ); 

  void resize(size_t size);

  iterator begin();
  iterator end();

  size_t size() const ; 
  bool get(size_t index) const; 
  void set(size_t index ) ; 
  void unset(size_t index);
  
  size_t count() const; 

  bitvector& operator|=(bitvector const& rhs );
  bitvector& operator&=(bitvector const& rhs );
  bitvector operator~() const ;

  bitvector symmetricDifference( bitvector const& rhs) const ; 

  bool operator<(bitvector const& rhs) const;
  bool operator>(bitvector const& rhs) const { return rhs < *this;  }

  bool operator<=(bitvector const& rhs) const { return (*this == rhs) ||  ( *this < rhs) ; }
  bool operator>=(bitvector const& rhs) const {return *this == rhs || *this > rhs; } 
  
  bool operator==(const bitvector &other) const ;
  bool operator!=(const bitvector &other) const { return !(*this == other); }

  friend std::ostream& operator<<(std::ostream& out, bitvector const& rhs );

  operator vector<size_t>() const ; 

public:                         // FRIENDS
  friend bitvector operator|( bitvector const& lhs, bitvector const& rhs);
  friend bitvector operator&( bitvector const& lhs, bitvector const& rhs);
  friend bitvector operator-(bitvector const& lhs, bitvector const rhs);

private:                        // METHODS
  static size_t elemSize(size_t size) ;
  void ensureCapacity(size_t c); 

private:
  vector<bv_elem> _array;
  size_t _usedSize;
};



// remove it again, the iterator sucks 
class bitvector::iterator
{
public:
  iterator(bitvector* rhs, size_t index);
  iterator& next() { ++_index; return *this; }
  iterator& operator++() { next(); return *this;  }

  bool get() const; 
  void set(bool val); 

private: 
  bitvector* _ref;
  size_t _index; 
} ;


inline bitvector::iterator::iterator(bitvector* rhs, size_t index)
  : _ref{rhs}
  , _index{index}
  {  }


inline bitvector::iterator bitvector::begin() { return  iterator(this, 0  );  }
inline bitvector::iterator bitvector::end() {return iterator(this, _array.size() ); }

#endif /* BITVECTOR_H */
