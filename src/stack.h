#ifndef _STACK_H
#define _STACK_H

#define INIT_STACK_LENGTH  10

#include "branch.h"

typedef struct _stack
{
  branch *content; 
  int length; 
  int index;   
} stack ; 

void createStack(stack **s); 
void clearStack(stack *s); 
void pushStack(stack *s, branch value); 
branch popStack(stack *s); 
int stackIsEmpty(stack *s); 
branch peekStack(stack *s); 
void printStack(stack *s ); 
void freeStack(stack **s); 
int stackLength(stack *s); 
#endif
