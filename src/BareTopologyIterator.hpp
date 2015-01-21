#ifndef BARETOPOLOGYITERATOR_H
#define BARETOPOLOGYITERATOR_H


class BareTopology::iterator
{
  friend class BareTopology;
  friend class Topology;

  using iterator_category = std::forward_iterator_tag ;
  using value_type = Link;
  using difference_type = size_t;
  using pointer = Link*;
  using reference = Link&;

  friend difference_type distance (iterator first, iterator last); 
    
public:
  iterator(BareTopology const*  ref, value_type l);
  iterator(  const iterator& rhs) = default;
  iterator& operator=(const iterator &rhs) = default;
  
  /** 
      @brief equals ->next and looses any traversal memory   
   */ 
  iterator neighbor() const;
  /** 
      @brief equals ->back and looses any traversal memory   
  */ 
  iterator opposite() const;
  /** 
      @brief advances the iterator by one element 
   */ 
  iterator& next();

  bool operator==(const iterator &other) const { return _ref == other._ref && _curLink == other._curLink; }
  bool operator!=(const iterator &other) const { return !(*this == other); }
  
  Link const& operator*() const {  return _curLink; }
  Link const *  operator->() const { return &_curLink ;  }

  iterator& operator++() { next(); return *this; }
  iterator operator+( difference_type n ) { auto cpy = *this; cpy.advance(n); return cpy; }
  /** 
      @brief advances the iterator by n elements 
   */ 
  iterator& advance(difference_type n);

  BareTopology const* getTopologyHandle() const { return _ref; }
  
  /** 
      @brief for DEBUG
   */ 
  Link get() const {return _curLink; }
  /** 
      @brief forgets the branches it already has traversed 
   */ 
  void reset() ;

  size_t getNumOuterNodesInSubtree() const ;
  size_t getNumLinksInSubtree() const;

  bool isValid() const;
  
  vector<node_id> findPathTo( iterator end) const ;


  static int ctr; 
  
private:
  vector<node_id> findPathTo_helper( iterator end) const;
  vector<node_id> findPathTo_helperWithBipartitions( bitvector const& endBip) const;

  bitvector getBipOrDummy() const ; 
private:
  BareTopology const* _ref; 
  value_type _curLink;
  std::stack<value_type> _descent;
  Link _first;
  bool _haveBothSides;
  size_t _numTraversed;       
};

#endif /* BARETOPOLOGYITERATOR_H */
