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


/**
   @brief compares if a pair of ints equals another pair. 

   @return TRUE on equality 
 */ 
boolean comparePair(int a1, int a2, int b1, int b2 ) 
{
  return (a1 == b1 && a2 == b2 ) 
    || (a1 == b2 && a2 == b1 ); 
}
