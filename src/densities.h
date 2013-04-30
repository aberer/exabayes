/**
   @file densities.h
   
   @brief contains densities and distributions
*/ 


#ifndef _DENSITIES_H
#define _DENSITIES_H

/** @brief density for dirichlet distribution with parameters "alphas" at point "values" */ 
double densityDirichlet(double *values, double *alphas, int length); 


/** @brief the gamma function */ 
double gammaFunction(double alpha); 


/** @brief the beta function */ 
double betaFunction(double *alpha, int length); 

#endif
