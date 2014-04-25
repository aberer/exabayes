#include "CommRequestImpl.hpp"

#include "system/extensions.hpp"

CommRequest::Impl::Impl(std::vector<char> array)
{
}


CommRequest::Impl::~Impl()
{
} 


bool CommRequest::Impl::isServed()
{
  return true; 
} 


std::vector<char> CommRequest::Impl::getArray() const
{
  return _array; 
} 
