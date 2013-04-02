#include "misc-utils.h" 

void swpDouble(double *a, double *b)
{
  double tmp = *b; 
  *b = *a; 
  *a = tmp; 
}

void swpInt(int *a, int *b)
{
  int tmp = *b;
  *b = *a; 
  *a = tmp; 
}

