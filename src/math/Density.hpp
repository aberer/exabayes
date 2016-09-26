/**
 * @file densities.h
 *
 * @brief contains densities and distributions
 */


#ifndef _DENSITIES_H
#define _DENSITIES_H

#include "extensions.hpp"

#include <vector>


// ____________________________________________________________________________
double                        exponentialDensity(
    double value,
    double lambda);


/** @brief density for dirichlet distribution with parameters "alphas" at point
 * "values" */
/* double densityDirichlet(double *values, double *alphas, int length);  */


/** @brief the gamma function */
// double gammaFunction(double alpha);

namespace Density {
// ____________________________________________________________________________
log_double                    lnDirichlet(
    std::vector<double>       values,
    const std::vector<double>&alphas);

// ____________________________________________________________________________
log_double                    lnExponential(
    double value,
    double lambda);

// ____________________________________________________________________________
log_double                    lnGamma(
    double x,
    double alpha,
    double beta);

// ____________________________________________________________________________
log_double                    lnWeibull(
    double x,
    double lambda,
    double k);

// ____________________________________________________________________________
double                        gammaFunction(
    double alpha);
}
#endif
