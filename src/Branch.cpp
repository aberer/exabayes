#include "Branch.hpp"
#include "parameters/AbstractParameter.hpp"
#include <limits>


template<>
void Branch<std::vector<double>>::deserialize( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   

  assert(0); 
  // TODO lengths 
}

template<>
void Branch<double>::deserialize( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   
  length =  cRead<decltype(length)>(in);
} 


template<>
void Branch<void>::deserialize( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   
}

template<>
void Branch<double>::serialize( std::ostream &out)  const
{
  cWrite<decltype(thisNode)>(out, thisNode); 
  cWrite<decltype(thatNode)>(out, thatNode); 
  cWrite<decltype(length)>(out, length); 
}

template<>
void Branch<std::vector<double>>::serialize( std::ostream &out)  const
{
  assert(0); 
}


template<>
void Branch<void>::serialize( std::ostream &out)  const
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
