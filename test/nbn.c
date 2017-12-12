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
    tile_n = 16;
    tile_m = 16;
    int tiled = 2*tile_m;
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
    
    for (i=0; i<tile_n && i<tile_m; i++)
    {
        remain -= 1;
        omp_set_nested(1);
        omp_set_num_threads((int)ceil(remain/(double)2));
        
        printf("loop %d , remain:%d , threads num:%d\n",i,remain,(int)(ceil(remain/2)));

        #pragma omp parallel shared(A,remain)
        {
            int thread_num,thread_id;
            thread_id = omp_get_thread_num();
            thread_num = omp_get_num_threads();

            mat *H;
            mat *Q;

            get_QR(A->m[thread_id*2*tile_m + i+i*tile_m] , &H, &Q);

            omp_set_num_threads(remain);
            #pragma omp parallel shared(A,H)
            {
                int sub_num,sub_id;
                sub_id = omp_get_thread_num();
                sub_num = omp_get_num_threads();
               #pragma omp critical
               {
                   printf("this is thread %d, update %d\n",sub_id,thread_id*2*tile_m+i+i*tile_m+sub_id);
               }
                mat * tmp_matrix = matrix_mul2(H , A->m[thread_id*2*tile_m+i+i*tile_m+sub_id]);
                matrix_free(A->m[thread_id*2*tile_m+i+sub_id]);
                A->m[thread_id*2*tile_m+i+sub_id] = tmp_matrix;
            }

            #pragma omp barrier
            matrix_free(H);

            if (i+(thread_id+1)*2 <= tile_n)
            {
                mat* merged = tile_merge(A,thread_id*2*tile_m+i+i*tile_m,thread_id*2*tile_m+i+i*tile_m+tile_m);
                get_QR(merged, &H, &Q);
                matrix_free(merged);
                omp_set_num_threads(tile_m - i);
                
                #pragma omp parallel shared(A,H,i)
                {
                    int sub_num,sub_id;
                    sub_id = omp_get_thread_num();
                    sub_num = omp_get_num_threads();

                    #pragma omp single
                    {
                        printf("this is thread %d mul on two %d,%d\n",sub_id,thread_id*tiled+i+i*tile_m+sub_id, thread_id*tiled+i+i*tile_m+sub_id+tile_m);
                    }
                    mat* new_merged = tile_merge(A, thread_id*tiled+i+i*tile_m+sub_id, thread_id*tiled+i+i*tile_m+sub_id+tile_m);
                    mat* tmp_matrix = matrix_mul2(H, new_merged);

                    matrix_free(new_merged);
                    int index[2] = {thread_id*tiled+i+i*tile_m+sub_id, thread_id*tiled+i+i*tile_m+sub_id+tile_m};

                    set_mul_value(A,tmp_matrix,index,2);

                    matrix_free(tmp_matrix);

                }
                matrix_free(H);
            }
            #pragma omp barrier
        }
        if (remain>2)
        {
           int tmp_count = 1;
            for (j=(int)ceil(remain/4.0);j>0;j=j/2)
            {
                tmp_count *= 2;
                printf("enter the j loop %d\n",j);
                omp_set_num_threads(j);
                #pragma omp parallel shared (A,j)
                {
                    int thread_id,thread_num;
                    thread_id = omp_get_thread_num();
                    thread_num = omp_get_num_threads();
                    int up = i*tile_m+i+thread_id*tile_m*2*tmp_count;
                    int down = i*tile_m+i+thread_id*tile_m*2*tmp_count+tile_m*tmp_count;
                    mat *H,*Q;
                    #pragma omp single
                    {
                        printf("this is %d merge on %d,%d\n",thread_id,up,down);
                    }
                    if (down<tile_m*tile_n)
                    {
                        mat* merged = tile_merge(A,up,down);
                        get_QR(merged, &H, &Q);
                        matrix_free(merged);
                    }
                        omp_set_num_threads(remain);

                        #pragma omp parallel shared(A,H)
                        {
                            if (down<tile_m*tile_n)
                            {
                                int sub_id, sub_num;
                                sub_id = omp_get_thread_num();
                                sub_num = omp_get_num_threads();
                                mat* new_merged = tile_merge(A,up+sub_id,down+sub_id);
                                mat* tmp_matrix = matrix_mul2(H, new_merged);
                                matrix_free(new_merged);
                                int index[2] = {up+sub_id,down+sub_id};
                           // #pragma omp critical
                            //{
                               // printf("This is thread %d, mul on two %d,%d\n",thread_id,up+sub_id,down+sub_id); 
                                set_mul_value(A,tmp_matrix,index,2);

                                matrix_free(tmp_matrix);
                            }
                        }
                        if (down<tile_m*tile_n)
                            matrix_free(H);
                }
            }
        } 
    }
    show_tile_matrix(A);
    return 1;
}
