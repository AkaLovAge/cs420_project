#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tile_matrix.h"
#include "matrix.h"

//initialize the tile matrix
//row: row number of original matrix 
//col: col number of original matrix 
//sub_row divide row into sub_row 
//sub_col divide col into sub_col
//the matrix is divide into:
// ...1...2...3...4...5...
// ...6...7...8...9...10..

mat_tile* tile_init(int row, int col, int sub_row, int sub_col)
{
    if (row%sub_row || col % sub_col)
    {
        printf("please enter row and col that can be divided by number of sub-chunk\n");
        exit(-1);
    }
    int n = row / sub_row;
    int m = col / sub_col;
    int size = sub_row * sub_col;

    mat_tile * new_ma = malloc(sizeof(mat_tile));
    new_ma->m = malloc(sizeof(mat*)*size);

    int i,j;
    for (i=0;i<size;i++)
    {
        new_ma->m[i] = matrix_init(n,m);

    }
    new_ma->row = sub_row;
    new_ma->col = sub_col;
    new_ma->sub_row = n;
    new_ma->sub_col = m;
    new_ma->num = sub_row * sub_col;

    return new_ma;
}

// free the tile matrix 

void matrix_tile_free(mat_tile* m)
{
    int i,j,k;
    for (i=0;i<m->num;i++)
    {
        matrix_free(m->m[i]);
    }
    free(m->m);
    free(m);
}

// set multiple matrix value to a tile_matrix , the matrix should be vertical

void set_mul_value(mat_tile* mt, mat* m, int index[], int num)
{
    int i,j,k;
    int row = mt->sub_row;
    int col = mt->sub_col;

    for (i=0;i<num;i++)
    {
        for (j=0;j<row;j++)
        {
            for(k=0;k<col;k++)
            {
                (mt->m[index[i]])->m[j][k] = m->m[i*row+j][k];
            }
        }
    }
}

// show tile matrix 

void show_tile_matrix(mat_tile* m)
{
    int i;
    for (i=0;i<m->num;i++)
    {
        printf("the num %d chunk:\n",i);
       
        matrix_show(m->m[i]);
    }
}

// convert regular matrix to tiled one

void matrix2tiled(mat_tile* m, double** rm )
{
    int i,j,k;
    int tmp1,tmp2;

    for (i=0;i<m->num;i++)
    {
        tmp1 = (i/m->row)*m->sub_row;
        tmp2 = (i%m->col)*m->sub_col;
        for (j=0; j<m->sub_row;j++)
        {
            for(k=0; k<m->sub_col; k++)
            {
                (m->m[i])->m[j][k] = rm[tmp1+j][tmp2+k];
            }
        }
    }
}

// copy a tile_matrix 

mat_tile* tile_matrix_copy(mat_tile* m)
{
    mat_tile* new_m = tile_init((m->row)*(m->sub_row),(m->col)*(m->sub_col), m->sub_row, m->sub_col);
        
    int sub_num = (m->sub_row)*(m->sub_col);
    int i,j,k;

    for (i=0;i<sub_num;i++)
    {
        for(j=0;j<m->sub_row;j++)
        {
            for(k=0;k<m->sub_col;k++)
            {
                (new_m->m[i])->m[j][k] = (m->m[i])->m[j][k];
            }
        }
    }
    return new_m;
}

//get nth colum from group matrix 
/*
double* get_col(mat* a, int n, int num1, int num2)
{
    int row = m->sub_row;
    double *tmp = malloc(sizeof(double)*row*2);
    int i;

    for (i=0;i<row;i++)
    {
        tmp[i] = (a->m[num1])->m[i][n];
        tmp[i+row] = (a->m[num2])->m[i][n];
    }

    return tmp;
}*/
mat* tile_merge(mat_tile* m, int num1, int num2)
{
    int row = m->sub_row;
    int col = m->sub_col;

    mat* new_ma = matrix_init(row*2,col);
    mat* mat1 = m->m[num1];
    mat* mat2 = m->m[num2];
    
    int i,j;

    for (i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            new_ma->m[i][j] = mat1->m[i][j];
            new_ma->m[i+row][j] = mat2->m[i][j];
        }
    }
    return new_ma;
}

// get H matrix for two matrix, the two matrix must be verticle
/*mat* tile_getH (mat_tile *m, int num1, int num2)
{
    if (num2 != num2 + m->col)
    {
        printf("the two matrix input is not verticle\n");
        exit(-1);
    }else{
        int i,j,k;
        int m_size = (m->sub_row > m->sub_col) ? m->sub_col:m->sub->row;
        mat* q[m_size];
        int row = m->sub_row;
        int col = m->sub_col;
        
        mat* tmp11 = m->m[num1];
        mat* tmp21 = m->m[num2];
        mat* tmp12,tmp22;

        double *e = malloc(sizeof(double)*row*2);

        double *n_col1;
        double n_norm;

        for (i=0;i<m_size;i++)
        {
            tmp12 = matrix_reflecor(tmp11,i);
            if (tmp12 != m->m[num1]) matrix_free(tmp11);
            tmp11 = tmp12;

            n_col = get_col2(m,i,num1,num2);
            n_norm = vnorm(n_col,row*2);
            n_norm = (n_col[i]<0) ? -n_norm:n_norm;

            for (j=0;j<row*2;j++) e[j]=0;
            e[i] = 1;

            vector_mul(e, n_norm, row*2);
            vector_add(n_col,e,row*2);

            n_norm = vnorm(n_col, row*2);
            vector_mul(n_col, (double)1/n_norm, row*2);

            q[i] = I_mul(n_col, row);

            tmp12 = matrix_mul2(q[i],tmp1)
        }
    }
}*/
