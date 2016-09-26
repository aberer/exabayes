#ifndef ARITHMETICS_HPP
#define ARITHMETICS_HPP

#include <vector>
#include "common.h"

namespace Arithmetics {
// ____________________________________________________________________________
/**
 * @brief gets the mean of the data
 */
double                    getMean(
    const std::vector<double>&data);

// ____________________________________________________________________________
/**
 *  @brief gets the variance of the data
 */
double                    getVariance(
    const std::vector<double>&data);

// ____________________________________________________________________________
/**
 *  @brief gets the potential scale reduction factor of data (data per chain)
 */
double                    PRSF(
    const std::vector<std::vector<double> >data);

// ____________________________________________________________________________
/**
 * @brief gets the n-th percentile of the data
 */
double                    getPercentile(
    double             percentile,
    std::vector<double>data);

// ____________________________________________________________________________
/**
 *  @brief gets the effective sampling size of the data
 */
double                    getEffectiveSamplingSize(
    const std::vector<double>& data);

// ____________________________________________________________________________
/**
 *  @brief gets the pearson correlation coefficient between two samples
 */
double                    getPearsonCorrelationCoefficient(
    const std::vector<double>&sampleA,
    const std::vector<double>&sampleB);

// ____________________________________________________________________________
double                    getAutoCorrelation(
    const std::vector<double>&data,
    nat                       lag);

// ____________________________________________________________________________
double                    getCoefficientOfVariation(
    const std::vector<double>&data);

// ____________________________________________________________________________
/**
 *  @brief calculates the sample skewness of data
 */
double                    getSkewness(
    const std::vector<double>&data);

// ____________________________________________________________________________
double                    getKahansSum2(
    const std::vector<double>&x);

// ____________________________________________________________________________
double                    getKahanSum(
    const std::vector<double>&x);
}
#endif
