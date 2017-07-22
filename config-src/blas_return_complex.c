
#if defined( MKL_ILP64 ) || defined( ILP64 )
    typedef long long myint;
#else
    typedef int myint;
#endif

#if defined( LOWERCASE )
    #define FORTRAN_NAME( lower, UPPER ) lower
#elif defined( UPPERCASE )       
    #define FORTRAN_NAME( lower, UPPER ) UPPER
#else                            
    #define FORTRAN_NAME( lower, UPPER ) lower ## _
#endif

#ifdef __cplusplus
    #define EXTERN_C extern "C"
#else
    #define EXTERN_C
#endif

#include <stdio.h>
#include <complex.h>

#define zdotc FORTRAN_NAME( zdotc, ZDOTC )

EXTERN_C
double _Complex zdotc( const myint* n,
                       const double _Complex* x, const myint* incx,
                       const double _Complex* y, const myint* incy );

int main( int argc, char** argv )
{
    myint n = 5, ione = 1;
    double _Complex x[5] = { 1, 2, 3, 4, 5 };
    double _Complex y[5] = { 5, 4, 3, 2, 1 };
    double _Complex expect = 35;
    double _Complex result = zdotc( &n, x, &ione, y, &ione );
    myint okay = (result == expect);
    printf( "zdotc result %.2f, expect %.2f, %s\n",
            creal(result), creal(expect), (okay ? "ok" : "failed"));
    return ! okay;
}
