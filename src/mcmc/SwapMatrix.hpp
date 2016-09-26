#ifndef _SWAP_MATRIX_H
#define _SWAP_MATRIX_H

#include "common.h"
#include "SuccessCounter.hpp"

#include "Serializable.hpp"

#include <vector>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//                                SWAP MATRIX                                //
///////////////////////////////////////////////////////////////////////////////
class SwapMatrix : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    explicit SwapMatrix(
        size_t numChains);
    // ________________________________________________________________________
    SwapMatrix(
        SwapMatrix&& rhs) = default;
    // ________________________________________________________________________
    SwapMatrix&                                    operator=(
        const SwapMatrix&rhs)  = default;
    // ________________________________________________________________________
    SwapMatrix&                                    operator=(
        SwapMatrix&&rhs)  = default;
    // ________________________________________________________________________
    SwapMatrix(
        const SwapMatrix& rhs)  = default;
    // ________________________________________________________________________
    void                                           update(
        nat  a,
        nat  b,
        bool acc);
    // ________________________________________________________________________
    const SuccessCounter&                          getCounter(
        nat a,
        nat b) const;
    // ________________________________________________________________________
    virtual void                                   deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                   serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    std::vector<SuccessCounter>                    getMatrix() const
    {return matrix; }
    // ________________________________________________________________________
    SwapMatrix                                     operator+(
        const SwapMatrix& rhs) const;
    // ________________________________________________________________________
    SwapMatrix&                                    operator+=(
        const SwapMatrix& rhs);
    // ________________________________________________________________________
    friend void                                    swap(
        SwapMatrix& a,
        SwapMatrix& b);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    size_t                                         mapToIndex(
        nat a,
        nat b) const;
    // ________________________________________________________________________
    friend  std::ostream&                          operator<<(
        std::ostream&     out,
        const SwapMatrix& rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<SuccessCounter> matrix;
    size_t                      numEntries;
};


#endif
