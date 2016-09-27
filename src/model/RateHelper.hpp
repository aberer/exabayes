/**
 *  @file RateHelper.hpp
 *
 *  @brief various functions helping to deal with rates and substition matrices
 *
 *  @notice using a class for this may appear over-engineered. But
 *  conversion between rate formats is not unproblematic, probably a
 *  Kahan summation may come in handy at some point.
 *
 */


#ifndef ABSTRACT_RATE_PROPOSER
#define ABSTRACT_RATE_PROPOSER

#include "common.h"

#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
//                                RATE HELPER                                //
///////////////////////////////////////////////////////////////////////////////
class RateHelper
{
public:
    RateHelper(){}
    // ________________________________________________________________________
    static void                                   convertRelativeToLast(
        std::vector<double>&values);
    // ________________________________________________________________________
    static double                                 convertToSum1(
        std::vector<double>&values);
    // ________________________________________________________________________
    static void                                   convertToGivenSum(
        std::vector<double>&values,
        double              givenSum);
    // ________________________________________________________________________
    static std::vector<double>                    getScaledValues(
        std::vector<double>values,
        double             scParameter);
    // ________________________________________________________________________
    static std::vector<nat>                       extractIndices(
        nat                       num,
        nat                       numRates,
        const std::vector<double>&rates);
    // ________________________________________________________________________
    static void                                   insertRates(
        nat                 num,
        nat                 numRates,
        std::vector<double>&rates,
        std::vector<double>&partRates);
    // ________________________________________________________________________
    static std::vector<double>                    extractSomeRates(
        nat                 num,
        nat                 numRates,
        std::vector<double>&rates);
    // ________________________________________________________________________
    static nat                                    numStateToNumInTriangleMatrix(
        int numStates);
};


#endif

