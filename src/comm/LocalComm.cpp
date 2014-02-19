#include "comm/LocalComm.hpp"
#include <cassert>


LocalComm::LocalComm(std::unordered_map<tid_t,int> tid2rank)
  : _messages(tid2rank.size())
  , _tid2LocCommIdx(tid2rank)
  , _colors(_tid2LocCommIdx.size(),0)
  , _ranks(_tid2LocCommIdx.size(), 0)
  , _size (_tid2LocCommIdx.size())
{
  nat ctr = 0; 
  for(auto &r : _ranks)
    r = ctr++; 
}

LocalComm::LocalComm(const LocalComm& rhs) 
  : _messages(rhs._messages)
  , _tid2LocCommIdx(rhs._tid2LocCommIdx)
  , _colors(rhs._colors)
  , _ranks(rhs._ranks)
  , _size(rhs._size)
  , _asyncMessages(rhs._asyncMessages)
{
}
  


LocalComm::LocalComm(LocalComm &&rhs) 
  : _messages(std::move(rhs._messages))
  , _tid2LocCommIdx(std::move(rhs._tid2LocCommIdx))
  , _colors(std::move(rhs._colors))
  , _ranks(std::move(rhs._ranks))
  , _size(std::move(rhs._size))
  , _asyncMessages(std::move(rhs._asyncMessages))
{
}

LocalComm& LocalComm::operator=(LocalComm rhs )
{
  swap(*this, rhs);
  return *this;
}


std::ostream& operator<<(std::ostream& out, const LocalComm& rhs)
{
  out << "{" << rhs.getColor() << "," << rhs.getRank() << "}"; 
  // for(nat i = 0; i < rhs._tid2LocCommIdx.size(); ++i)
  //   out << "{" <<  rhs._colors[i] << "," << rhs._ranks[i] << "}" ; 
  return out; 
}


int LocalComm::getNumThreads() const 
{
  return _ranks.size();
}

int LocalComm::size() const 
{
  return _size; 
}


int LocalComm::getRank() const
{ 
  return _ranks.at(_tid2LocCommIdx.at(MY_TID));  
} 


bool LocalComm::isValid() const  
{
  // TODO evil, but reasonable 
  bool result = true; 
  for(nat i = 0; i < _colors.size() ; ++i)
    for(nat j = i + 1 ; j < _colors.size() ; ++j)
      result &= (_ranks[i] != _ranks[j] ||  _colors[i] != _colors[j] ) ; 
  
  auto color2num =  std::unordered_map<int,int>{};
  for(nat i =0; i < _colors.size() ; ++i)
    ++color2num[_colors[i]] ; 
  int aNum = begin(color2num)->second; 
  for(auto &elem : color2num)
    assert(aNum == elem.second); 
  assert(_size == aNum); 

  return result; 
}


LocalComm LocalComm::split(const std::vector<int> &color, const std::vector<int> &rank)  const
{
  auto result = *this; 
  
  auto uniqCols = std::unordered_set<int>{}; 
  uniqCols.insert(begin(color), end(color)); 
  auto minCol = *(std::min_element(begin(uniqCols), end(uniqCols))); 
  auto tColors = color; 

  if(rank[0] > 0 )
    {
      auto tRanks = rank; 
      for(auto &r : tRanks)
	r -= rank[0]; 
      result.setRanks(tRanks); 
      result.setColors(tColors); 
    }
  else 
    {
      result.setRanks(rank);
      result.setColors(tColors); 
    }
  
  auto col2occ = std::unordered_map<int,int>{}; 
  for(auto c : result._colors)
    ++col2occ[c]; 
  result._size = begin(col2occ)->second; 

  result.isValid(); 
 
  return result; 
}

 
void swap(LocalComm& lhs, LocalComm& rhs )
{
  using std::swap; 
  swap(lhs._messages, rhs._messages); 
  swap(lhs._tid2LocCommIdx, rhs._tid2LocCommIdx); 
  swap(lhs._colors, rhs._colors); 
  swap(lhs._ranks, rhs._ranks); 
  swap(lhs._size, rhs._size); 
} 


int LocalComm::getIdx(int col, int rank) const 
{
  int result = -1; 
  for(int i = 0; i < getNumThreads(); ++i)
    {
      if(_colors[i] == col && _ranks[i] == rank)
	return i; 
    }
  assert(0); 
  return result; 
}



void LocalComm::waitAtBarrier()
{
  // TODO efficiency. But it is not used currently anyway...
  auto arr = std::vector<uint8_t>{1};
  arr = allReduce<uint8_t>(arr);
  assert(arr.size() == 1 && arr[0] == size()); 
}



static void dummyFun()
{
  
}

bool LocalComm::haveThreadSupport() const 
{
  std::cout << "TODO proper check for pthread support " << std::endl; 
  return true; 
}



bool LocalComm::checkAsyncMessage( int tag )  const 
{
  return _asyncMessages.checkMessage(tag); 
}
