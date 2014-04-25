#include "IncompleteMesh.hpp"
#include <cassert>
#include <algorithm>
#include <tuple>

IncompleteMesh::IncompleteMesh(nat size, nat runDimSize, nat chainDimSize)
  : _runDimSize (runDimSize)
  , _chainDimSize ( chainDimSize)
  , _globalSize( size) 
{
}


std::tuple<nat,nat> IncompleteMesh::getElementsPerDimension(nat total, nat dimSize) const 
{
  nat cartProcs = total  / dimSize; 
  nat remainderProcs = total % dimSize ; 
  nat numMoreProc = remainderProcs; 

  assert(numMoreProc * (cartProcs + 1 )  +  (  dimSize - numMoreProc ) * cartProcs == total); 
  return std::make_tuple(  cartProcs, dimSize - numMoreProc ); 
}


nat IncompleteMesh::getProcsInMyDim(nat rank, nat total, nat dimSize ) const 
{
  auto elemPerDim = getElementsPerDimension(total, dimSize); 
  auto procs = std::get<0>(elemPerDim); 
  auto numWithFewProcs = std::get<1>(elemPerDim); 
  return  ( rank < ( dimSize - numWithFewProcs ) * (procs+1)  )? procs + 1  : procs ; 
}


nat IncompleteMesh::getMyCoord(nat rank, nat total, nat dimSize ) const 
{
  auto elemPerDim = getElementsPerDimension(total, dimSize); 
  auto procs = std::get<0>(elemPerDim); 
  auto numWithFew  = std::get<1>(elemPerDim); 
  auto numWithMany = dimSize - numWithFew; 

  if( rank <  numWithMany  * ( procs + 1) )
    return rank / (procs + 1 ); 
  else 
    {
      rank -=  numWithMany * (procs + 1); 
      return  numWithMany  + rank /  procs; 
    }
}


nat IncompleteMesh::getRankInMyDim(nat rank, nat total, nat dimSize) const 
{
  nat procs = 0, numWithFew = 0; 
  std::tie(procs, numWithFew) = getElementsPerDimension(total, dimSize); 
  nat numWithMany = dimSize - numWithFew;  
    
  if(rank < numWithMany * (procs + 1 ))
    return rank % (procs+1); 
  else 
    {
      rank -=  numWithMany * (procs + 1); 
      return rank % procs; 
    }
}


std::array<nat,3> IncompleteMesh::getCoordinates(nat rank) const
{
  auto result = std::array<nat,3>{{0,0,0}}; 

  result[0] = getMyCoord(rank, _globalSize, _runDimSize); 
    
  auto numProcMyRun = getProcsInMyDim(rank, _globalSize, _runDimSize); 
  auto rankInRun = getRankInMyDim(rank, _globalSize, _runDimSize); 
  result[1] = getMyCoord(rankInRun, numProcMyRun, _chainDimSize); 

  auto numProcMyChain = getProcsInMyDim( rankInRun, numProcMyRun, _chainDimSize ); 
  auto rankInChain = getRankInMyDim(rankInRun, numProcMyRun, _chainDimSize ); 

  result[2] = getMyCoord(rankInChain, numProcMyChain, numProcMyChain); 
    
  return result; 
}


std::ostream& operator<<(std::ostream& out, const IncompleteMesh& rhs)
{
  for(nat i = 0; i < rhs._globalSize; ++i)
    {
      auto coord = rhs.getCoordinates(i); 
      out << "rank:\t" << i << "\t|\t"<< std::get<0>(coord) << "\t" << 
      	std::get<1>(coord) << "\t" << 
      	std::get<2>(coord) <<  std::endl; 
    }
  return out; 
}


nat IncompleteMesh::getRankFromCoordinates( std::array<nat,3> coords) const 
{
  auto result = 0; 
  auto resultRun = 0; 
  auto resultChain = 0; 

  assert(coords[0] < _runDimSize && coords[1] < _chainDimSize); 
  auto elemPerDim = getElementsPerDimension(_globalSize, _runDimSize); 

  // add ranks for run 
  auto procs = std::get<0>(elemPerDim); 
  auto numWithFew = std::get<1>(elemPerDim);
  auto numWithMany = _runDimSize - numWithFew; 
  if(coords[0] < numWithMany)
    resultRun += coords[0] * (procs + 1 ); 
  else 
    {
      resultRun += numWithMany * (procs + 1); 
      coords[0] -= numWithMany; 
      resultRun += coords[0] * procs; 
    }

  // add ranks
  auto numProcMyRun = getProcsInMyDim(resultRun, _globalSize, _runDimSize); 
  auto rankInRun = getRankInMyDim(resultRun, _globalSize, _runDimSize); 
  assert(rankInRun == 0); 
  elemPerDim = getElementsPerDimension(numProcMyRun, _chainDimSize); 
  procs = std::get<0>(elemPerDim); 
  numWithFew = std::get<1>(elemPerDim);
  numWithMany = _chainDimSize - numWithFew; 
    
  if(coords[1] < numWithMany)
    resultChain += coords[1] * (procs+1); 
  else 
    {
      resultChain += numWithMany * (procs+1); 
      coords[1] -= numWithMany; 
      resultChain += coords[1] *  procs; 
    }
    
  result += resultRun + resultChain + coords[2]; 

  return result; 
}


nat IncompleteMesh::getNumRanksInDim(nat runBatchId, nat chainBatchId) const 
{
  assert(runBatchId < _runDimSize); 
  assert(chainBatchId < _chainDimSize); 
  
  auto baseRank = getRankFromCoordinates( {{ runBatchId, chainBatchId, 0 }}); 
  nat otherRank = 0;  
  
  if(chainBatchId + 1 < _chainDimSize)
    {
      otherRank = getRankFromCoordinates( {{ runBatchId, chainBatchId + 1 , 0 }}); 
    }
  else if(runBatchId + 1 < _runDimSize)
    {
      otherRank = getRankFromCoordinates( {{ runBatchId + 1 , 0 , 0}} ); 
    }
  else 				
    {
      otherRank = _globalSize; 
    }

  assert(baseRank  < otherRank); 

  return otherRank - baseRank; 
}   
