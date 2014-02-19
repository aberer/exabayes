#ifndef _INCOMPLETE_MESH_HPP
#define _INCOMPLETE_MESH_HPP

#include "common.h"
#include <iostream>
#include <array>


class IncompleteMesh
{
public: 
  IncompleteMesh(nat size, nat runDimSize, nat chainDimSize); 
  /** 
      @brief gets the rank from coordinates in this mesh 
   */ 
  nat getRankFromCoordinates( std::array<nat,3> coords) const ; 
  /** 
      @brief gets the coordinates of a process in this mesh from the rank 
   */ 
  std::array<nat,3> getCoordinates(nat rank) const; 
  /** 
      @brief gets the effective size of a dimension 
      
      in other words the total number of ranks assigned to something 
   */ 
  nat getNumRanksInDim(nat runBatchId, nat chainBatchId) const ;   

  nat getRunDimSize() const {return _runDimSize; }
  nat getChainDimSize() const {return _chainDimSize; }
  

private: 			// METHODS 
  /** 
      @brief gets number of elements per dimension and number of
      batches for which this is applicable

      Notice that dimSize - result[1] will have result[0] + 1 elements
      
      @param total -- total  number of elements 
      @paarm dimSize -- size of the dimension
   */ 
  std::tuple<nat,nat> getElementsPerDimension(nat total, nat dimSize) const ; 
  nat getProcsInMyDim(nat rank, nat total, nat dimSize ) const ; 
  nat getMyCoord(nat rank, nat total, nat dimSize ) const ; 
  nat getRankInMyDim(nat rank, nat total, nat dimSize) const ; 
  friend std::ostream& operator<<(std::ostream& out, const IncompleteMesh& rhs); 

private: 			// ATTRIBUTES
  nat _runDimSize; 
  nat _chainDimSize; 
  nat _globalSize;
};  

#endif
