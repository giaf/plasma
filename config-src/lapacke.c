
#include <stdio.h>

#ifdef HAVE_MKL
    #include <mkl_lapacke.h>
#else
    #include <lapacke.h>
#endif

int main( int argc, char** argv )
{
    int i, n = 2, info = 0;
    double A[2*2] = { 16, 4,   -1, 5 };
    double L[2*2] = {  4, 1,   -1, 2 };
    double work[1];
    info = LAPACKE_dpotrf_work( LAPACK_COL_MAJOR, 'L', n, A, n );
    if (info != 0) {
        printf( "dpotrf failed: info %d\n", info );
        return 1;
    }
    for (i = 0; i < n*n; ++i) {
        if (A[i] != L[i]) {
            printf( "dpotrf failed: A[%d] %.2f != L[%d] %.2f\n",
                    i, A[i], i, L[i] );
            return 1;
        }
    }
    printf( "dpotrf ok\n" );
    return 0;
}
