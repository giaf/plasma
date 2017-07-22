#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "plasma.h"

int main()
	{

	printf("\nciao init\n");

	int ii;

	// problem description

	plasma_enum_t uplo = PlasmaLower;

	int m = 1024;

	int Am, An;

	Am = m;
	An = m;

	int lda = Am;

	double *A = (double*)malloc((size_t)lda*An*sizeof(double));
	for(ii=0; ii<m*m; ii++)
		A[ii] = 0.0;
	for(ii=0; ii<m; ii++)
		A[ii*(lda+1)] = 1.0;

	omp_set_num_threads(4);

	plasma_init();

	// plasma tuning parameters
	plasma_set(PlasmaTuning, PlasmaDisabled);
	plasma_set(PlasmaNb, 64);

	int rep;
	int nrep = 10;
	struct timeval tv0, tv1;

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		plasma_dpotrf(uplo, m, A, lda);
		}

	gettimeofday(&tv1, NULL); // stop

	double time = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);
	double Gflops = 1.0/3.0*m*m*m*1e-9/time;

	printf("\ntime   = %e\n", time);
	printf("\nGflops = %f\n", Gflops);

	plasma_finalize();

	free(A);

	printf("\nciao end\n");

	return 0;

	}


