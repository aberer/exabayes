#include <assert.h>
#include "axml.h"

#include "stack.h" 

void createStack(stack **s)
{
  *s = (stack*)exa_calloc(1,sizeof(stack)); 
  (*s)->content = (branch*)exa_calloc(INIT_STACK_LENGTH, sizeof(branch)); 
  (*s)->length = INIT_STACK_LENGTH; 
  (*s)->index = 0;  
}


void clearStack(stack *s)
{
  s->index = 0;   
}


void pushStack(stack *s, branch value)
{
  if(s->index == s->length)
    {
      s->content = (branch*)exa_realloc(s->content,s->length * 2  * sizeof(branch));
      s->length *= 2;
    }
  
  s->content[s->index++] = value; 
}

branch popStack(stack *s)
{
  if(s->index == 0)
    {
      assert(0); 
      branch b = {0,0}; 
      return b; 
    }
  else 
    {
      s->index--; 
      return s->content[s->index]; 
    }
}


int stackIsEmpty(stack *s)
{
  return s->index == 0; 
}


int stackLength(stack *s)
{
  return s->index; 
}



branch peekStack(stack *s)
{
  if(s->length == 0)
    {
      assert(0); 
      branch b = {0,0};
      return b; 
    }
  else 
    return s->content[s->index-1]; 
}


void printStack(stack *s )
{
  printf("content: ");
  for(int i =0; i < s->index; ++i)
    printf("{%d,%d},", s->content[i].thisNode, s->content[i].thatNode); 
  printf("\n"); 
}


void freeStack(stack **s)
{
  exa_free((*s)->content); 
  exa_free(*s); 
  *s = NULL; 
}
