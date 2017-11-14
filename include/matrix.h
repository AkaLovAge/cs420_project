#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int row,col;
    double ** m;
} mat;

mat* matrix_init(int row, int col);

mat* matrix_copy(int row, int col,mat* input);

double* get_col(mat* input,int n);

void matrix_transpose(mat * input);

void vector_add(double *a, double *b,int n);

void vector_mul(double *a, double s, int n);

void matrix_add(mat* a, mat* b);

mat* matrix_mul(mat* input1, mat* input2);

mat * matrix_reflector(mat *a, int n);

double vnorm(double* input,int num);

mat *I_mul(double *v,int n);

void matrix_free(mat* input);

void matrix_show(mat *m);


