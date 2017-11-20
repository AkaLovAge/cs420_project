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
        tmp1 = (i%m->row)*m->sub_row;
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


