#ifndef _INCOMPLETE_MESH_HPP
#define _INCOMPLETE_MESH_HPP

#include "common.h"

#include <iostream>
#include <array>

///////////////////////////////////////////////////////////////////////////////
//                              INCOMPLETE MESH                              //
///////////////////////////////////////////////////////////////////////////////
class IncompleteMesh
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    IncompleteMesh(
        size_t size,
        size_t runDimSize,
        size_t chainDimSize);
    // ________________________________________________________________________
    /**
     *  @brief gets the rank from coordinates in this mesh
     */
    nat                                           getRankFromCoordinates(
        std::array<nat, 3>coords) const;
    // ________________________________________________________________________
    /**
     *  @brief gets the coordinates of a process in this mesh from the rank
     */
    std::array<nat, 3>                            getCoordinates(
        nat rank) const;
    // ________________________________________________________________________
    /**
     *  @brief gets the effective size of a dimension
     *
     *  in other words the total number of ranks assigned to something
     */
    nat                                           getNumRanksInDim(
        nat runBatchId,
        nat chainBatchId) const;
    // ________________________________________________________________________
    size_t                                        getRunDimSize() const
    {return _runDimSize; }
    // ________________________________________________________________________
    size_t                                        getChainDimSize() const
    {return _chainDimSize; }

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    // METHODS
    /**
     *  @brief gets number of elements per dimension and number of
     *  batches for which this is applicable
     *
     *  Notice that dimSize - result[1] will have result[0] + 1 elements
     *
     *  @param total -- total  number of elements
     *  @paarm dimSize -- size of the dimension
     */
    std::tuple<size_t, size_t>                    getElementsPerDimension(
        size_t total,
        size_t dimSize) const;
    // ________________________________________________________________________
    size_t                                        getProcsInMyDim(
        nat    rank,
        size_t total,
        size_t dimSize) const;
    // ________________________________________________________________________
    nat                                           getMyCoord(
        nat    rank,
        size_t total,
        size_t dimSize) const;
    // ________________________________________________________________________
    nat                                           getRankInMyDim(
        nat    rank,
        size_t total,
        size_t dimSize) const;
    // ________________________________________________________________________
    friend std::ostream&                          operator<<(
        std::ostream&         out,
        const IncompleteMesh& rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    size_t _runDimSize;
    size_t _chainDimSize;
    size_t _globalSize;
};


#endif
