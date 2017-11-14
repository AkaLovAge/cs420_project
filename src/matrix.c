#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <matrix.h>
// initialize the mat*

mat* matrix_init(int row, int col)
{
    mat* new_ma = malloc(sizeof(mat));
    new_ma->m = malloc(sizeof(double*)*row);
    int i;
    for (i=0; i<row; i++)
    {
        new_ma->m[i] = calloc(sizeof(double),col);
    }
    new_ma->row = row;
    new_ma->col = col;
    return new_ma;
}

// copy input mat to new mat

mat* matrix_copy(int row, int col,mat* input)
{
    if(input->row != row || input->col != col)
    {
        printf("input matrix is not %d * %d\n",row,col);
        exit(-1);
    }

    mat* new_ma;
    new_ma = matrix_init(row,col);
    int i,j;

    for (i=0; i<row; i++)
    {
        for (j=0; j<col; j++)
        {
            new_ma->m[i][j] = input->m[i][j];
        }
    }

    return new_ma;
}

// get nth column of matrix, return a malloc double list

double* get_col(mat* input,int n)
{
    int row = input->row;
    double *tmp = malloc(sizeof(double)*row);
    int i;

    for (i=0; i<row; i++)
    {
        tmp[i] = input->m[i][n];
    }
    
    return tmp;
}

// mat m[i][j] = m[j][i]

void matrix_transpose(mat * input)
{
    int i,j;
    if (input->row == input->col)
    {
        for (i=0;i<input->row;i++)
        {
            for(j=0;j<i;j++)
            {
                double tmp = input->m[i][j];
                input->m[i][j] = input->m[j][i];
                input->m[j][i] = tmp;
            }
        }
    }else{
        printf("can not transpose matrix\n");
    }
}

// vector a = a + b
   
void vector_add(double *a, double *b,int n)
{
    if(sizeof(a) != sizeof(b))
    {
        printf("two vector can not add up\n");
    }else{
        int i;
        for (i=0;i<n;i++)
        {
            a[i] += b[i];
        }
    }
}

// vector a = a * s

void vector_mul(double *a, double s, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
        a[i]=a[i] * s;
    }
}

// matrix a = a + b

void matrix_add(mat* a, mat* b)
{
    if (a->row != b->row || a->col != b-> col)
    {
        printf("two matrix can not add up.\n");
    }else{
        int i,j;

        for (i=0;i<a->row;i++)
        {   
            for(j=0;j<a->col;j++)
            {                   
                a->m[i][j] = a->m[i][j]+b->m[i][j];
            }                                   
        }                                               
    }                                                       
}                   

//mat a * mat b

mat* matrix_mul(mat* input1, mat* input2)
{
    int i,j,k;
    mat* result = matrix_init(input1->row,input2->col);
    if (input1->col == input2->row)
    {
        for (i=0;i<input1->row;i++)
        {
            for(j=0;j<input2->col;j++)
            {
                double tmp=0;
                for(k=0;k<input1->col;k++)
                {
                    tmp = tmp+(input1->m[i][k])*(input2->m[k][j]);
                }
                result->m[i][j] = tmp;
            }
        }
    }else{
        printf("col of m1 and row of m2 is not equal.\n");
    }

    return result;
}

// generate [I,0][0,x]

mat * matrix_reflector(mat *a, int n)
{
    mat *result = matrix_init(a->row,a->col);
    int i,j,k;
    for (i=0;i<n;i++)
    {
        result->m[i][i] = 1;
    }
    for (j=n;j<a->row;j++)
    {
        for(k=n;k<a->col;k++)
        {
            result->m[j][k] = a->m[j][k];
        }
    }
    
    return result;
}

//||vector|| return double 

double vnorm(double* input,int num)
{
    int i;
    double tmp=0;
    for(i=0;i<num;i++)
    {
        tmp += input[i]*input[i];
    }

    return sqrt(tmp);
}

// I - 2*v*v^T

mat *I_mul(double *v,int n)
{
    mat *result = matrix_init(n,n);
    int i,j;
    for (i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            result->m[i][j] = -2*v[i]*v[j];
        }
     }

     for (i=0;i<n;i++)
        result->m[i][i] = 1+result->m[i][i];

    return result;
}

// free mat malloc 

void matrix_free(mat* input)
{
    if(input->row && input->col)
    {
        int i;
        for(i=0; i<input->row; i++)
        {
            free(input->m[i]);
        }
        free(input->m);
        free(input);
    }else{
        printf("can not free matrix\n");
        exit(-1);
    }
}

// show the mat 

void matrix_show(mat *m)
{
    int i,j;
    for (i=0;i<m->row;i++)
    {
        for(j=0;j<m->col;j++)
        {
            printf("%f,",m->m[i][j]);
        }
        printf("\n");
    }
}

