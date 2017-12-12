#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <matrix.h>

// initialize the mat*

mat* matrix_init(int row, int col)
{
    mat* new_ma = malloc(sizeof(mat));
    new_ma->m = malloc(sizeof(double*)*row);
    new_ma->m[0] = calloc(sizeof(double),row*col);
    int i;
    for (i=1; i<row; i++)
    {
        new_ma->m[i] = new_ma->m[0] + i*col;
    }
    new_ma->row = row;
    new_ma->col = col;
    return new_ma;
}

// set matrix a value to matrix b 

void set_value(mat *a, mat *b)
{
    if (a->row != b->row || a->col != b->col)
    {
        int i;
        for (i=0;i<(a->row)*(a->col);i++)
        {
            a->m[0][i] = b->m[0][i];
        }
    }
}
// copy input mat to new mat

mat* matrix_copy(int row, int col,mat* input)
{
    if(input->row != row || input->col != col)
    {
        printf("input matrix is not %d * %d\n",row,col);
        exit(-1);
    }

    mat* new_ma =  matrix_init(row,col);
    int i,j;

    for (i=0; i<row*col; i++)
    {
            new_ma->m[0][i] = input->m[0][i];
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
// mat a = mat a * mat b
void matrix_mul(mat *a, mat *b)
{
    int i,j,k;
    mat* result = matrix_init(a->row,b->col);
    if (a->col == b->row)
    {
         for (i=0;i<a->row;i++)
         {
             for(j=0;j<b->col;j++)
             {
                 double tmp = 0;
                 for(k=0;k<b->col;k++)
                 {
                     tmp = tmp + (a->m[i][k])*(b->m[k][j]);
                 }
                 result->m[i][j] = tmp;
             }
        }
    }
}

// regular mat a * mat b
mat* matrix_mul2(mat* input1, mat* input2)
{
    int i,j,k;
    mat* result = matrix_init(input1->row,input2->col);
  
    if (input1->col == input2->row)
    {
        for (i=0;i<input1->row;i++)
        {
            for(j=0;j<input2->col;j++)
            {
                double tmp = 0;
                for(k=0;k<input1->col;k++)
                {
                    tmp = tmp + (input1->m[i][k])*(input2->m[k][j]);
                }
                result->m[i][j] = tmp;
            }
         }
     }else{
         printf("col of m1 and row of m2 is not equal.\n"); 
     }
     return result;
}

//mat a * mat b, input cut b from col m1 to col m2
mat* matrix_mul4(mat* input1, mat* input2, int n1, int n2, int m1, int m2)
{
    if (input1->col == (n2 - n1+1))
    {    
        int i,j,k;
        mat* result = matrix_copy(input2->row, input2->col,input2);

        for (i=0;i<input1->row;i++)
        {
            for(j=m1;j<m2+1;j++)
            {
                double tmp=0;
                for(k=0; k<input1->col; k++)
                {
                    tmp = tmp+ (input1->m[i][k])*(input2->m[k+n1][j]);
                }
                result->m[i][j] = tmp;
            }
        }
        return result;
    }else{
        printf("col of m1 and row of m2 is not equal.\n");
    }
    return NULL;
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
        free(input->m[0]);
        free(input->m);
        free(input);

    }else{
        printf("can not free matrix\n");
        exit(-1);
    }
}

// only free the data in matrix 

//void matrix_data_free(mat *input)
//{
//    if (input-> && input ->)
//}

// householder return R and Q matrix, must input which block is included, from n1,m1 to n2,m2 

void get_QR(mat *m, mat **H, mat **Q)
{
    int i,j,k;
    int row = m->row;
    int col = m->col;
    int l_r = row , l_c = col;
    int p_size = 0;
    if (row == col)
    {
        p_size = (m->row < m->col)? (m->row-1): (m->col-1);
 //       l_r = row-1;
 //       l_c = col-1;

    }else{
        p_size = (row < col)? row:col;
    }
    
    mat *q[p_size];
    mat *tmp1 = m;
    mat *tmp2;

    double *e = malloc(sizeof(double)*row);
    double *n_col;
    double n_norm;

    for(i=0;i<l_r && i<l_c; i++)
    {
        tmp2 = matrix_reflector(tmp1,i);
        if (tmp1 != m) matrix_free(tmp1);

        tmp1 = tmp2;
        n_col = get_col(tmp1 , i);
 //        for (j=0;j<row;j++)
 //                        printf("%f,",n_col[j]);

   //                              printf("\n");
   //     for (j=0;j<row;j++)
   //         printf("%f,",n_col[j]);
   //     printf("\n");
        n_norm = vnorm(n_col,row);
        n_norm = (n_col[i]<0) ? -n_norm : n_norm;

        for(j=0; j<row; j++) {e[j]=0;}
        e[i] = 1;
        vector_mul(e,n_norm, row);
        vector_add(n_col,e,row);
        n_norm = vnorm(n_col,row);
        vector_mul(n_col, (double)1/n_norm, row);

        q[i] = I_mul(n_col, row);
      //  matrix_show(q[i]);   
        tmp2 = matrix_mul2(q[i],tmp1);
        free(n_col);
    
//        printf ("%d\n",i);
 //       matrix_show(q[i]);

        if (tmp1 != m) matrix_free(tmp1);

        tmp1 = tmp2;
     }
     free(e);
     matrix_free(tmp1);
     mat* tmp5 = q[0];
     mat* tmp4;
   //  matrix_show(tmp5);
     for(i=1; i<l_r && i<l_c; i++)
     {
        tmp4 = tmp5;
        tmp5=matrix_mul2(q[i], tmp5);
        matrix_free(tmp4);
     }
  //   matrix_show(tmp5);
     *H = matrix_copy(tmp5->row, tmp5->col,tmp5);
     matrix_free(tmp5);
}

// show the mat 

void matrix_show(mat *m)
{
    int i,j;
    printf("row:%d,col:%d\n",m->row,m->col);
    for (i=0;i<m->row;i++)
    {
        for(j=0;j<m->col;j++)
        {
            printf("%.20f,",m->m[i][j]);
        }
        printf("\n");
    }
}

