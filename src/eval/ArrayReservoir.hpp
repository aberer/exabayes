#ifndef REVERVE_ARRAYS_HPP
#define REVERVE_ARRAYS_HPP

#include "common.h"

#include <list>
#include <memory>
#include <map>
#include <unordered_map>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
//                              ARRAY RESERVOIR                              //
///////////////////////////////////////////////////////////////////////////////
class ArrayReservoir
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ArrayReservoir(
        bool freeOlds);
    // ________________________________________________________________________
    ~ArrayReservoir();
    // ________________________________________________________________________
    ArrayReservoir(
        const ArrayReservoir& rhs) = delete;
    // ________________________________________________________________________
    ArrayReservoir&                         operator=(
        const ArrayReservoir& rhs) = delete;
    // ________________________________________________________________________
    double*                                 allocate(
        size_t requiredLength);
    // ________________________________________________________________________
    void                                    deallocate(
        double* array);

#ifdef _DEVEL
    std::tuple<uint64_t,
               uint64_t>                    getUsedAndUnusedBytes()
    const;
#endif

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::unordered_map<double*, size_t>   _usedArrays;
    std::map<size_t, std::list<double*> > _unusedArrays;

    static const double                   numGammaCats;
    static const double                   thresholdForNewSEVArray;
    bool                                  _freeOlds;
};


#endif
