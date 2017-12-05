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
    tile_n = 4;
    tile_m = 4;
    int i,j,k;

    double **a = malloc(sizeof(double*)*N);
    for (i=0;i<N;i++)
    {
        a[i] = malloc(sizeof(double)*M);
        for (j=0;j<M;j++)
            a[i][j] = (double)rand()/((double) RAND_MAX/100);
    }

    mat_tile* A = tile_init(N,M,tile_n,tile_m);
    matrix2tiled(A,a);
    
    int remain = tile_n+1;
    int thread1,thread2;
    
   // for (i=0; i<tile_n && i<tile_m; i++)
   // {
    i=0;
        remain -= 1;
        omp_set_nested(1);
        omp_set_num_threads((int)(ceil(remain/2)));
        
        #pragma omp parallel shared(A)
        {
            int thread_num,thread_id;
            thread_id = omp_get_thread_num();
            thread_num = omp_get_num_threads();

            mat *H;
            mat *Q;

            get_QR(A->m[thread_id*8 + i+i*tile_m] , &H, &Q);

            omp_set_num_threads(tile_m - i);
            #pragma omp parallel shared(A,H)
            {
                int sub_num,sub_id;
                sub_id = omp_get_thread_num();
                sub_num = omp_get_num_threads();
                
                mat * tmp_matrix = matrix_mul2(H , A->m[thread_id*8+i+sub_id]);
                matrix_free(A->m[thread_id*8+i+sub_id]);
                A->m[thread_id*8+i+sub_id] = tmp_matrix;
            }

            #pragma omp barrier
            matrix_free(H);

            if (i+(thread_id+1)*2 <= tile_n)
            {
                mat* merged = tile_merge(A,thread_id*8+i+i*tile_m,thread_id*8+i+i*tile_m+4);
                get_QR(merged, &H, &Q);
                matrix_free(merged);
                omp_set_num_threads(tile_m - i);
                
                #pragma omp parallel shared(A,H)
                {
                    int sub_num,sub_id;
                    sub_id = omp_get_thread_num();
                    sub_num = omp_get_num_threads();

                    mat* new_merged = tile_merge(A, thread_id*8+i+i*tile_m+sub_id, thread_id*8+i+i*tile_m+sub_id+4);
                    mat* tmp_matrix = matrix_mul2(H, new_merged);

                    matrix_free(new_merged);
                    int index[2] = {thread_id*8+i+i*tile_m+sub_id, thread_id*8+i+i*tile_m+sub_id+4};

                    set_mul_value(A,tmp_matrix,index,2);

                    matrix_free(tmp_matrix);

                }
            }

        }
        show_tile_matrix(A);
    //}
}
