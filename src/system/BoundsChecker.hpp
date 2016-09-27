#ifndef _BOUNDS_CHECKER
#define _BOUNDS_CHECKER

#include "BranchLength.hpp"
#include "BranchLengths.hpp"

#include <vector>

///////////////////////////////////////////////////////////////////////////////
//                               BOUNDS CHECKER                              //
///////////////////////////////////////////////////////////////////////////////
class BoundsChecker
{
public:
    // ________________________________________________________________________
    static const double zMin, zMax,
                        rateMax, rateMin,
                        alphaMin, alphaMax,
                        freqMin;
    // ________________________________________________________________________
    static bool                    checkFrequencies(
        const std::vector<double>&freqs);
    // ________________________________________________________________________
    static bool                    checkBranch(
        const BranchLength&branch);
    // ________________________________________________________________________
    static bool                    checkBranch(
        const BranchLengths&branch);
    // ________________________________________________________________________
    static bool                    checkRevmat(
        const std::vector<double>&rates);
    // ________________________________________________________________________
    static bool                    checkAlpha(
        double alpha);
    // ________________________________________________________________________
    static void                    correctAlpha(
        double&alpha);
    // ________________________________________________________________________
    static void                    correctBranch(
        BranchLengths&branch);
    // ________________________________________________________________________
    static void                    correctBranch(
        BranchLength&branch);
    // ________________________________________________________________________
    static void                    correctRevMat(
        std::vector<double>&rates);
    // ________________________________________________________________________
    static void                    correctFrequencies(
        std::vector<double>&frequencies);
}


;


#endif
