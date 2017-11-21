#ifndef TILE_MATRIX_H
#define TILE_MATRIX_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"

typedef struct mat_tile{
    int row, col;
    int sub_row, sub_col;
    int num;
    mat **m;
} mat_tile;

mat_tile* tile_init(int row, int col, int sub_row, int sub_col);

void matrix_tile_free(mat_tile* m);

void show_tile_matrix(mat_tile* m);

void matrix2tiled(mat_tile* m, double **rm);

mat_tile* tile_matrix_copy(mat_tile* m);

mat* tile_merge(mat_tile* m, int num1, int num2);

void set_mul_value(mat_tile* mt, mat* m, int index[], int num);

#endif
