#ifndef ___RANDOMNESS_H
#define ___RANDOMNESS_H

#include "GlobalVariables.hpp"
#include "Serializable.hpp"

#include <iostream>
#include <vector>

/* TODO is this correct? */
#ifndef UINT64_C
#define UINT64_C
#endif

/*
 * we could mess around here a lot with 64-bit or the even cooler
 * alternative to threefry (aesni). For the time being that's hardly
 * worth it, threefry is already much better than the default RNG.
 */

#include "threefry.h"
#include "uniform.hpp"
#define exa_rand(c, k) threefry2x32(c, k)

using randKey_t = threefry2x32_key_t;
using randCtr_t = threefry2x32_ctr_t;

class TreeAln;

///////////////////////////////////////////////////////////////////////////////
//                                 RANDOMNESS                                //
///////////////////////////////////////////////////////////////////////////////
class Randomness : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    Randomness(
        randCtr_t seed);
    // ________________________________________________________________________
    void                                   setKey(
        randKey_t key);
    // ________________________________________________________________________
    randCtr_t                              generateSeed();
    // ________________________________________________________________________
    /**
     *  @brief sets the generation for the RNG.
     */
    void                                   rebaseForGeneration(
        uint64_t generation);
    // ________________________________________________________________________
    uint64_t                               getGeneration() const;
    // ________________________________________________________________________
    /**
     *  @brief draw a floating point random  number in [0,1)
     */
    double                                 drawRandDouble01();
    // ________________________________________________________________________
    nat                                    operator()();
    // ________________________________________________________________________
    double                                 drawRandBiUnif(
        double x);
    // ________________________________________________________________________
    /**
     *  @brief draw from a exponential distribution with parameter
     *  lambda
     */
    double                                 drawRandExp(
        double lambda);
    // ________________________________________________________________________
    /**
     * @brief draw a random number that is distributed uniformly on the
     * exponential scale
     */
    double                                 drawMultiplier(
        double multiplier);
    // ________________________________________________________________________
    /**
     *  @brief draw a random number around param for a given window
     */
    double                                 drawFromSlidingWindow(
        double param,
        double window);
    // ________________________________________________________________________
    randCtr_t                              getKey() const {return key; }
    // ________________________________________________________________________
    /**
     *  @brief draws a permutation
     */
    void                                   drawPermutation(
        int* perm,
        int  n);
    // ________________________________________________________________________
    /**
     * @brief draw r according to distribution given by weights.
     *
     * NOTE sum of weights is not required to be 1.0
     */
    // int drawSampleProportionally( double *weights, int numWeight );
    /**
     *  @brief draws a random number form a dirichlet distribution
     */
    std::vector<double>                    drawRandDirichlet(
        const std::vector<double>&alphas);
    // ________________________________________________________________________
    /**
     *  @brief draw a random number from a gamma distribution with shape
     *  parameter alpha and scale parameter beta
     */
    double                                 drawRandGamma(
        double alpha,
        double beta);
    // ________________________________________________________________________
    /**
     *  @brief This function should be called if the alphas for the dirichlet
     * distribution are given
     */
    void                                   drawRandDirichlet(
        std::vector<double>&      results,
        const std::vector<double>&alphas);
    // ________________________________________________________________________
    /**
     *  @brief This function should be called if the expected values for the
     * dirichlet distribution are given
     */
    void                                   drawDirichletExpected(
        std::vector<double>&      results,
        const std::vector<double>&mean,
        double                    scale);
    // ________________________________________________________________________
    /**
     *  @brief draws integer uniformly from [0,n)
     */
    template<typename T>
    T                                      drawIntegerOpen(
        T upperBound);
    // ________________________________________________________________________
    double                                 drawRandWeibull(
        double lambda,
        double k);
    // ________________________________________________________________________
    /**
     *  @brief draw integer uniformly from [0,n]
     */
    template<typename T>
    T                                      drawIntegerClosed(
        T upperBound);
    // ________________________________________________________________________
    nat                                    drawGeometric(
        double prop);
    // ________________________________________________________________________
    nat                                    drawBinomial(
        double prop,
        nat    trials);
    // ________________________________________________________________________
    virtual void                           deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                           serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&    out,
        const Randomness&rhs);
    // ________________________________________________________________________
    nat                                    max()
    {return std::numeric_limits<nat>::max(); }
    // ________________________________________________________________________
    nat                                    min()
    {return std::numeric_limits<nat>::min(); }

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                                   incrementNoLimit();

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    randCtr_t ctr;
    randKey_t key;
};


#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wtautological-compare"


template<typename T>
T                                          Randomness::drawIntegerOpen(
    T upperBound)
{
    assert(int64_t(upperBound) > 0 && upperBound <=
           std::numeric_limits<nat>::max());
    return drawIntegerClosed(upperBound - 1);
}

template<typename T>
T                                          Randomness::drawIntegerClosed(
    T upperBound)
{
    assert(upperBound >= 0 && upperBound <= std::numeric_limits<nat>::max());

    if (upperBound == 0)
        return 0;
    else
    {
        randCtr_t r = exa_rand(key, ctr);
        ++ctr.v[0];
        return r.v[0] % (upperBound + 1);
    }
}

#endif
