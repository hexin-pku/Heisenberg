#include <stdio.h>

//lapacke headers
#include "lapacke.h"
#include "lapacke_config.h"
#include "lapacke_utils.h"

extern lapack_int LAPACKE_dgesv( int matrix_order, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );

int main()
{
    int N = 4;
    double A[16] = {  1,  2,  3,  1,
                      4,  2,  0,  2,
                     -2,  0, -1,  2,
                      3,  4,  2, -3};
    double B[8] = {  6,  2,  1,  8,
                     1,  2,  3,  4};
    int ipiv[4];
    int n = N;
    int nrhs = 2;
    int lda = N;
    int ldb = N;

    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR,n,nrhs,A,lda,ipiv,B,ldb);
    printf("info:%d\n",info);
    if(info==0)
    {
        int i = 0;
        int j = 0;
        for(j=0;j<nrhs;j++)
        {
            printf("x%d\n",j);
            for(i=0;i<N;i++)
                printf("%.6g \t",B[i+j*N]);
            printf("\n");
        }
    }
}
