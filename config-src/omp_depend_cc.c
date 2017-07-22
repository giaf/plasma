
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
    int last = 1000;
    for (int iter = 0; iter < 100; ++iter) {
        // inserts last/10 tasks that update x
        #pragma omp parallel
        {
            for (int i = 0; i <= last; i += 10) {
                #pragma omp task depend(inout:x[0:n])
                task( n, x, i );
            }
        }
        // verify that updates worked
        for (int i = 0; i < n; ++i) {
            int expect = last + i;
            if (x[i] != expect) {
                printf( "openmp task depend failed, x[%d] = %d, expected %d (iter %d)\n",
                        i, x[i], expect, iter );
                return 1;
            }
        }
    }
    printf( "openmp task depend seems ok\n" );
    return 0;
}
