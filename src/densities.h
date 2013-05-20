/**
   @file densities.h
   
   @brief contains densities and distributions
*/ 


#ifndef _DENSITIES_H
#define _DENSITIES_H

#include <vector>
using namespace std; 

double exponentialDensity(double value, double lambda); 


/** @brief density for dirichlet distribution with parameters "alphas" at point "values" */ 
double densityDirichlet(double *values, double *alphas, int length); 


/** @brief the gamma function */ 
double gammaFunction(double alpha); 


/** @brief the beta function */ 
double betaFunction(double *alpha, int length); 

double densityDirichletLog(double *values, double *alphas, int length); 
double densityDirichletWrapper(vector<double> values, vector<double> alphas, bool logScale); 
double exponentialDistribution(double value, double lambda); 

#endif
