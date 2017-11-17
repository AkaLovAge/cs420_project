#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "timers.h"
#include "mpi.h"
#include "omp.h"

void MPI_check(char str[], int error)
{
    if (error != MPI_SUCCESS )
    {
        printf("error when call %s, exit.\n",str);
        MPI_Finalize();
        exit(-1);
    }
}

int main (int argc, char* argv[])
{
    int i,j,k;
    mat* A = matrix_init(N,M);
    
    for (i=0;i<N;i++)
    {
        for (j=0; j<M; j++)
        {
            A->m[i][j] = (double)rand()/((double)RAND_MAX/10);
        }
    }

    #pragma omp parallel
    {
        int threads_num,thread_id;
        thread_id = omp_get_thread_num();
        threads_num = omp_get_num_threads();
     
     // getting the H of first 2d matrix 
        #pragma omp single 
        {

        }

    }

}
