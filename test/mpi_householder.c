#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "timers.h"
#include "mpi.h"
//#include "mpi_writer.h"

//double A[N][M];
void MPI_Check(char str[],int error)
{
    if (error != MPI_SUCCESS )
    {
        printf("error happend when call MPI %s, exit.\n",str);
        MPI_Finalize();
        exit(-1);
    }
}
int main(int argc, char* argv[])
{
    int n,m,m_size;
    int mpi_error;
    int size,rank;
    MPI_Status status;
    MPI_File fh;
    MPI_Offset filesize;
    MPI_Comm comm;

// MPI initial
   
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("The rank%d is on, total has %d, in size\n",rank,size);

// get the column and row for each processor and find the size of each of them 
    n = N/size;
    m = M;
    printf("The rank%d is on, total has %d, in size %d*%d\n",rank,size,n,m);

// initialize the matrix, needed to rewrite
    double tmp;
    int i,j,k;
    mat* A = matrix_init(n,m);
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
        {
            tmp = (float)rand()/((float)RAND_MAX/10);
            A->m[i][j] = tmp;
        }
    }

// go into householder
    mat *tmp1 = A;
    mat* tmp2;
    double *e = malloc(sizeof(double)*N);
    double n_norm;
    mat *q;
/*    for (i=0;i<n && i<m;i++)
    {

//each processor bcast and collect whole column 
        double col4e[N];
        double *tmp_col = get_col(A,i);
        mpi_error = MPI_Allgather(tmp_col,n,MPI_DOUBLE, col4e, n, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Check("allgather",mpi_error);

//calculate normalization         
        n_norm =vnorm (col4e,N);
        n_norm = (col4e[i]<0) ? -n_norm:n_norm;

// initailize e vector
        for(j =0; j<N;j++){e[j]=0;}
        e[i] = 1;
        vector_mul(e,n_norm,N);
        vector_add(col4e,e,N);
        
        n_norm = vnorm(col4e,N);
        vector_mul(col4e,(double)1/n_norm,N);
    
        q = I_mul(col4e,N);

        //tmp2 = matrix_mul(q,A);
        free(tmp_col);
        MPI_Barrier(MPI_COMM_WORLD);
   }*/
        MPI_Finalize();
    return 0;
}

