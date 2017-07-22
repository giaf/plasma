
#include <stdio.h>

#if __STDC_VERSION__ >= 199901L
    // supports C99
#else
    choke function();
#endif

int main( int argc, char** argv )
{
    printf( "hello, __STDC_VERSION__ = %ld\n", __STDC_VERSION__ );
    for (int i = 0; i < 10; ++i) {
        // pass
    }
    return 0;
}
