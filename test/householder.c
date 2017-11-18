#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <matrix.h>
#include "tile_matrix.h"
/*
void householder(mat *m, mat *R, mat *Q)
{
    int i,j,k;
    mat *q[m->row];
    mat *tmp1 = m;
    mat * tmp2;
    // some where need to be parallized
    double *e = malloc(sizeof(double)*(m->row));
    double *n_col;
    double n_norm;
    for (i=0;i<m->row && i<m->col; i++)
    {
    	
    	tmp2 = matrix_reflector(tmp1,i);

	if (tmp1 != m) matrix_free(tmp1);
	tmp1 = tmp2;
	
	n_col = get_col(tmp1,i);
	n_norm = vnorm(n_col,m->row);

	n_norm = (m->m[i][i]<0) ? -n_norm:n_norm;
	for(j =0; j<m->row;j++){e[j]=0;}
	e[i] = 1;
	
	vector_mul(e,n_norm,m->row);
	vector_add(n_col,e,m->row);

	n_norm = vnorm(n_col,m->row);
	vector_mul(n_col,(double)1/n_norm,m->row);
	
	q[i]=I_mul(n_col,m->row);
	tmp2 = matrix_mul(q[i],tmp1);
	free(n_col);
	if (tmp1!=m) matrix_free(tmp1);
	tmp1 = tmp2;
    }
    free (e);
    matrix_free(tmp1);
    mat *tmp3=m;
    mat *tmp4;
    mat* tmp5=q[0];
    for (i=0;i<m->row && i<m->col;i++)
    {
    	if (i!=0) matrix_mul(tmp5,q[i]);
    	tmp4 = matrix_mul(q[i],tmp3);
	if (tmp3!=m) matrix_free(tmp3);
	tmp3 = tmp4;
	matrix_show(tmp3);
    }
    R = tmp3;
    matrix_show(tmp5);
}
*/
int main(int argc, char* argv[])
{
    double in[8][8] = {
    	{ -1, -1, 1,23,45,3,4,4},
	    {  1, 3, 3,3,4,5,2,5},
	    { -1, -1, 5,2,3,4,5,6},
	    { 1, 3, 7,2,3,5,6,7},
        { -1, -1, 1,23,45,3,4,4},
        {  1, 3, 3,3,4,5,2,5},
        { -1, -1, 5,2,3,4,5,6},
        { 1, 3, 7,2,3,5,6,7},
    };
    double **in1 = malloc(sizeof(double *)*8);
    int i,j,k;
    for (i=0;i<8;i++)
    {
        in1[i] = malloc(sizeof(double)*8);
        for(j=0;j<8;j++)
        {
            in1[i][j] = in[i][j];
        }
    }
    mat_tile * m = tile_init(8,8,2,2);

    matrix2tiled(m, in1);

    show_tile_matrix(m);

    mat_tile* m1 = tile_matrix_copy(m);

    show_tile_matrix(m1);
    matrix_tile_free(m);
    matrix_tile_free(m1);
}

