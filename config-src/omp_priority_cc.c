
#include <stdio.h>

void task( int n, int* x, int id )
{
    for (int i = 0; i < n; ++i) {
        x[i] = id + i;
    }
}

int main( int argc, char** argv )
{
    int n = 1000, x[1000] = { 0 };
    #pragma omp parallel
    {
        #pragma omp task depend(inout:x[0:n]) priority(1)
        task( n, x, 0 );

        #pragma omp task depend(inout:x[0:n]) priority(2)
        task( n, x, 100 );
    }
    for (int i = 0; i < n; ++i) {
        if (x[i] != 100 + i) {
            printf( "openmp task priority failed, x[%d] = %d, expected %d\n",
                    i, x[i], 100 + i );
            return 1;
        }
    }
    printf( "openmp task priority ok\n" );
    return 0;
}
