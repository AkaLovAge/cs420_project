#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "tile_matrix.h"
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
        printf("this is thread %d in total %d\n",thread_id, threads_num); 
     // getting the H of first 2d matrix 
     
        mat *H;
        mat *Q;
        #pragma omp single 
        {
            int m = M/threads_num -1;
            int n = N/threads_num -1;
            printf("the %d thread go for %d,%d\n",thread_id,n,m);

            get_QR_mn(A, 0, n, 0, m, &H, &Q);
        }
    //    #pragma omp for schedule (auto) private(i,j,k)
    //    {
//            for(i =0; i<M; i++)
  //              for(j =0; j<n;j++)
    //                for(k=)
    //    }
        printf("this is %d\n",thread_id);
        matrix_show(H);

    }

    return 1;
}
