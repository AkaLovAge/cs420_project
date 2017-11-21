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
    int tile_n,tile_m;
    
    tile_n = 2;
    tile_m = 2;

    int i,j,k;
    double ** a = malloc(sizeof(double*)*N);
    for (i=0;i<N;i++)
    {
        a[i] = malloc(sizeof(double)*M);
        for (j=0;j<M;j++)
            a[i][j] = (float)rand()/((float)RAND_MAX/10);
    }

    mat_tile* A = tile_init(N,M,tile_n,tile_m);
    matrix2tiled(A,a);
    
    mat *H;
    mat *Q;
   
    for (i=0; i<tile_n && i<tile_m; i++)
    {    
        get_QR(A->m[i+i*tile_m], &H, &Q);

        mat *tt = matrix_mul2(H,A->m[i+i*tile_m]);
   
        omp_set_num_threads(tile_m-i);
        #pragma omp parallel shared(A,H)
        {
            int threads_num,thread_id;
            thread_id = omp_get_thread_num();
            threads_num = omp_get_num_threads();
        
/*        #pragma omp critical
        {
            printf ("this is thread %d\n",thread_id);
            matrix_show(A->m[thread_id]);
        }
*/
            mat* tmp_matrix = matrix_mul2(H,A->m[i+thread_id]);
            matrix_free(A->m[i+thread_id]);
            A->m[i+thread_id] = tmp_matrix;
        }
        matrix_free(H);

        if(i != tile_n-1 || i != tile_m-1)
        {

    // get the H on verticle 
    //        omp_set_nested(1);
    
            mat* merged = tile_merge(A,0,0+A->row);
    
            get_QR(merged,&H,&Q);
            matrix_free(merged);

            #pragma omp parallel shared(A,H)
            {
                int threads_num,thread_id;
                thread_id = omp_get_thread_num();
                threads_num = omp_get_num_threads();
        
                mat* new_merged = tile_merge(A,thread_id,thread_id+A->row);
                mat* tmp_matrix = matrix_mul2(H,new_merged);
                matrix_free(new_merged);
                int index[2] = {thread_id,thread_id+A->col};
        
                set_mul_value(A,tmp_matrix,index,2);
        
                matrix_free(tmp_matrix);
            }
        }
    }
    show_tile_matrix(A);
//#pragma omp parallel shared(A,H)
//    {
//        int num1,id1;
//        id1 = omp_get_thread_num();
//        num1 = omp_get_num_threads();
/*
        #pragma omp barrier 

        #pragma omp single 
        {
            printf("this is thread %d\n",thread_id);
            show_tile_matrix(A);
        }
*/
      //  #pragma omp parallel num_threads(2)
      //  {
      //      int num2,id2;
      //      id2 = omp_get_thread_num();
      //     num2 = omp_get_num_threads();
      //  }
//    }
    return 1;
}
