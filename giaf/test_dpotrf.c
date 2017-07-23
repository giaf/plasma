#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "plasma.h"

int main()
	{

	int ii;

	// problem description

	plasma_enum_t uplo = PlasmaLower;
	char cl = 'l';
	int info = 0;

	int nb = 64;

	int mmax = 4096;

	omp_set_num_threads(4);

//	openblas_set_num_threads(1);

	plasma_init();

	plasma_set(PlasmaTuning, PlasmaDisabled);

	printf("\nmm\tnb\ttime\tGflops\n");

	int mm;
	for(mm=nb; mm<=mmax; mm+=nb)
		{

		int m = mm;

		int Am, An;

		Am = m;
		An = m;

		int lda = Am;

		double *A = (double*)malloc((size_t)lda*An*sizeof(double));
		for(ii=0; ii<m*m; ii++)
			A[ii] = 1e-3;
		for(ii=0; ii<m; ii++)
			A[ii*(lda+1)] = 2.0;
		for(ii=0; ii<m-1; ii++)
			A[1+ii*(lda+1)] = 1.0;

		// plasma tuning parameters
		plasma_set(PlasmaNb, nb);

		struct timeval tvi, tv0, tv1;
		double tvtmp, tvmin = 1e12;

		gettimeofday(&tvi, NULL); // start

		int rep;
		int nrep = 10;
		while(1)
			{

			for(ii=0; ii<m*m; ii++)
				A[ii] = 1e-3;
			for(ii=0; ii<m; ii++)
				A[ii*(lda+1)] = 2.0;
			for(ii=0; ii<m-1; ii++)
				A[1+ii*(lda+1)] = 1.0;

			gettimeofday(&tv0, NULL); // start

			plasma_dpotrf(uplo, m, A, lda);
//			dpotrf_(&cl, &m, A, &lda, &info);

			gettimeofday(&tv1, NULL); // stop

			tvtmp = (tv1.tv_sec-tv0.tv_sec)+(tv1.tv_usec-tv0.tv_usec)*1e-6;
			tvmin = tvtmp<tvmin ? tvtmp : tvmin;

			tvtmp = (tv1.tv_sec-tvi.tv_sec)+(tv1.tv_usec-tvi.tv_usec)*1e-6;
			
			if(tvtmp>1)
				break;

			}

	//	d_print_mat(m, m, A, lda);

		double time = tvmin;
		double Gflops = 1.0/3.0*m*m*m*1e-9/time;

		printf("%d\t%d\t%e\t%f\n", mm, nb, time, Gflops);

		free(A);

		}

	plasma_finalize();

	return 0;

	}


