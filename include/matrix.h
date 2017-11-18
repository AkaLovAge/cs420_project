#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int row,col;
    double ** m;
} mat;

mat* matrix_init(int row, int col);

mat* matrix_copy(int row, int col,mat* input);

double* get_col(mat* input,int n, int n1, int n2);

void matrix_transpose(mat * input);

void vector_add(double *a, double *b,int n);

void vector_mul(double *a, double s, int n);

void matrix_add(mat* a, mat* b);

mat* matrix_mul4(mat* input1, mat* input2,int n1, int n2, int m1, int m2);

mat* matrix_mul2(mat* input1, mat* input2);

mat * matrix_reflector(mat *a, int n);

double vnorm(double* input,int num);

mat *I_mul(double *v,int n);

void matrix_free(mat* input);

void matrix_show(mat *m);

void get_QR_mn(mat*m, int n1, int n2, int m1, int m2, mat **R, mat **Q);

#endif
