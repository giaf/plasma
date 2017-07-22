
#include <stdio.h>
#include <omp.h>

int main( int argc, char** argv )
{
    int x[10], nt;

    #pragma omp parallel
    nt = omp_get_num_threads();

    #pragma omp parallel for
    for (int i = 0; i < 10; ++i) {
        x[i] = i;
    }
    printf( "openmp x[0]=%d, nt=%d\n", x[0], nt );
    return 0;
}
