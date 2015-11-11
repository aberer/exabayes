/**
   @file densities.h
   
   @brief contains densities and distributions
*/ 


#ifndef _DENSITIES_H
#define _DENSITIES_H

#include <vector>

double exponentialDensity(double value, double lambda); 


/** @brief density for dirichlet distribution with parameters "alphas" at point "values" */ 
/* double densityDirichlet(double *values, double *alphas, int length);  */


/** @brief the gamma function */ 
double gammaFunction(double alpha); 


double densityDirichletLog(std::vector<double> values, const std::vector<double> &alphas); 
double densityDirichlet(std::vector<double> values, const std::vector<double> &alphas); 

double exponentialDistribution(double value, double lambda); 

double logGammaDensity (double x, double alpha, double beta ); 
#endif
