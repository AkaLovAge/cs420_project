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
    double ** a = malloc(sizeof(double*)*N);
    for (i=0;i<N;i++)
    {
        a[i] = malloc(sizeof(double)*M);
        for (j=0;j<M;j++)
            a[i][j] = i;
    }
    mat_tile* A = tile_init(N,M,2,8);
    matrix2tiled(A,a);
    show_tile_matrix(A);
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
            //get_QR_mn(A, 0, n, 0, m, &H, &Q);
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
