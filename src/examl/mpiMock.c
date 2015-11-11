#include "mpiMock.h"
#include <string.h>
#include <assert.h>

#ifdef _SEQUENTIAL

void MPI_Abort(int tmp, int tmp2) {} 
void MPI_Comm_size(int tmp, int *mtp2) { *mtp2 = 1 ; }
void MPI_Init(int *argc, char ***argv) {} 
void MPI_Comm_rank(int tmp, int *tmp2){ *tmp2 = 0; }
void MPI_Finalize(){}
void MPI_Comm_split(int tmp , int tmp2, int tmp3 , int *tmp4){}
void MPI_Barrier(int tmp){}
void MPI_Gatherv(void *x, int y, MPI_Datatype z, void *a, void* b, void *c, MPI_Datatype d, int x1, MPI_Comm foo){}
void MPI_Bcast(void *x, int y, MPI_Datatype z, int a, MPI_Comm b){} 
void MPI_Reduce(void *x, void *y, int z, MPI_Datatype t, MPI_Operator o, int x3, MPI_Comm comm) 
{
  if(t == MPI_DOUBLE)
    memcpy((double*)y,(double*)x, z * sizeof(double)); 
  else if(t == MPI_INT)
    memcpy((int*)y,(int*)x,z * sizeof(int)); 
  else 
    assert(0); 

}
void MPI_Gather(void *x, int y, MPI_Datatype z, void *a,  int someSize, MPI_Datatype d, int x1, MPI_Comm foo){}
void MPI_Allreduce(void *send, void *recv, int tmp , MPI_Datatype type, MPI_Operator op, MPI_Comm comm)
{
  MPI_Reduce(send, recv,tmp, type, op, 0,comm);
} 


#endif
