#include "Branch.hpp"
#include "parameters/AbstractParameter.hpp"
#include <limits>

BranchEqualFlag operator|( BranchEqualFlag a, BranchEqualFlag b) 
{
  return static_cast<BranchEqualFlag>(static_cast<int>(a) | static_cast<int>(b)); 
}


BranchEqualFlag operator&( BranchEqualFlag a, BranchEqualFlag b)
{
  return static_cast<BranchEqualFlag>(static_cast<int>(a) & static_cast<int>(b)); 
}


// std::ostream& operator<<(std::ostream &out, const Branch<void>& br)
// {
//   out << "(" << br.getPrimNode() << "/" << br.getSecNode() << ")"; 
//   return out; 
// }

// std::ostream& operator<<(std::ostream &out, const Branch<double>& br)
// {
//   out << "(" << br.getPrimNode() << "/" << br.getSecNode() << ")"; 
//   out << ":"  << br.getLength(); 
//   return out; 
// }

// std::ostream& operator<<(std::ostream &out, const Branch<std::vector<double>>& br)
// {
//   out << "(" << br.getPrimNode() << "/" << br.getSecNode() << ")"; 
//   out << ":[" ; 
//   bool isFirst = true; 
//   for(auto &v : br.getLengths())
//     {
//       if(isFirst)
// 	isFirst = false; 
//       else 
// 	out << ","; 
//       out << v ; 
//     }
//   out << "]"; 
//   return out; 
// }






template<>
void Branch<std::vector<double>>::readFromCheckpoint( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   

  assert(0); 
  // TODO lengths 
}

template<>
void Branch<double>::readFromCheckpoint( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   
  length =  cRead<decltype(length)>(in);
} 


template<>
void Branch<void>::readFromCheckpoint( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   
}

template<>
void Branch<double>::writeToCheckpoint( std::ostream &out)  const
{
  cWrite<decltype(thisNode)>(out, thisNode); 
  cWrite<decltype(thatNode)>(out, thatNode); 
  cWrite<decltype(length)>(out, length); 
}

template<>
void Branch<std::vector<double>>::writeToCheckpoint( std::ostream &out)  const
{
  assert(0); 
}


template<>
void Branch<void>::writeToCheckpoint( std::ostream &out)  const
{
  cWrite<decltype(thisNode)>(out, thisNode); 
  cWrite<decltype(thatNode)>(out, thatNode); 
}  



template<>
Branch<void> Branch<void>::getInverted() const
{
  auto result = Branch<void>(thatNode,thisNode); 
  return result; 
}

template<>
Branch<double> Branch<double>::getInverted() const
{
  auto result = Branch<double>(thatNode,thisNode); 
  result.length = this->length; 
  
  // if(std::is_same<BT, std::vector<double>>::value)
  //   result.lengths = this->lengths; 

  return result; 
} 

template<>
Branch<std::vector<double>> Branch<std::vector<double>>::getInverted() const
{
  auto result = Branch<std::vector<double>>(thatNode,thisNode); 
  result.lengths = this->lengths; 
  return result; 
}



template<>
Branch<double> Branch<std::vector<double>>::toOneLength(const AbstractParameter* param) const
{
  auto result = Branch<double>(thisNode,thatNode); 
  result.setLength(getLength(param)); 
  return result; 
} 
