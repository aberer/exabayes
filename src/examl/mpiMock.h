#ifndef _MPIMOCK
#define _MPIMOCK

/* mocks the mpi interface */

typedef int MPI_Comm ; 
typedef  int MPI_Datatype; 

#define MPI_INT 0 
#define MPI_DOUBLE 1

#define MPI_COMM_WORLD 0

typedef int MPI_Operator; 
#define MPI_SUM 0  


void MPI_Abort(int tmp, int tmp2) ; 
void  MPI_Comm_size(int tmp, int *mtp2) ;
void MPI_Init(int *argc, char ***argv) ; 
void MPI_Comm_rank(int tmp, int *tmp2);
void MPI_Finalize();
void MPI_Comm_split(int tmp , int tmp2, int tmp3 , int *tmp4);
void MPI_Barrier(int tmp);
void MPI_Gatherv(void *x, int y, MPI_Datatype z, void *a, void* b, void *c, MPI_Datatype d, int x1, MPI_Comm foo);
void MPI_Bcast(void *x, int y, MPI_Datatype z, int a, MPI_Comm b); 
void MPI_Reduce(void *x, void *y, int z, MPI_Datatype t, MPI_Operator o, int x3, MPI_Comm comm) ;
void MPI_Gather(void *x, int y, MPI_Datatype z, void *a,  int someSize, MPI_Datatype d, int x1, MPI_Comm foo);
void MPI_Allreduce(void *send, void *recv, int tmp , MPI_Datatype type, MPI_Operator op, MPI_Comm comm); 

#endif
