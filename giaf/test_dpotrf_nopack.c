#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "plasma.h"
#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_tuning.h"
#include "plasma_types.h"
#include "plasma_workspace.h"


int main()
	{

	int ii;

	// problem description

	plasma_enum_t uplo = PlasmaLower;
	char cl = 'l';
	int info = 0;

	int nb = 64;

	int mmax = 4096;
//	int mmax = 128;

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

		// Create tile matrix.
		plasma_desc_t dA;
		int retval;
		retval = plasma_desc_triangular_create(PlasmaRealDouble, uplo, nb, nb, m, m, 0, 0, m, m, &dA);

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

			plasma_context_t *plasma = plasma_context_self();

			// Set tiling parameters.
			int nb = plasma->nb;

			// Initialize sequence.
			plasma_sequence_t sequence;
			retval = plasma_sequence_init(&sequence);

			// Initialize request.
			plasma_request_t request;
			retval = plasma_request_init(&request);

//			// asynchronous block
			#pragma omp parallel
			#pragma omp master
				{
//				// Translate to tile layout.
				plasma_omp_dtr2desc(A, lda, dA, &sequence, &request);
				}
			// implicit synchronization

			gettimeofday(&tv0, NULL); // start

//			plasma_dpotrf(uplo, m, A, lda);

			// asynchronous block
			#pragma omp parallel
			#pragma omp master
				{
				// Call the tile async function.
				plasma_omp_dpotrf(uplo, dA, &sequence, &request);
				}
			// implicit synchronization

			gettimeofday(&tv1, NULL); // stop

			// asynchronous block
			#pragma omp parallel
			#pragma omp master
				{
				// Translate back to LAPACK layout.
				plasma_omp_ddesc2tr(dA, A, lda, &sequence, &request);
				}
			// implicit synchronization

			tvtmp = (tv1.tv_sec-tv0.tv_sec)+(tv1.tv_usec-tv0.tv_usec)*1e-6;
			tvmin = tvtmp<tvmin ? tvtmp : tvmin;

			tvtmp = (tv1.tv_sec-tvi.tv_sec)+(tv1.tv_usec-tvi.tv_usec)*1e-6;
			
			if(tvtmp>1)
				break;

			}

		// Free matrix A in tile layout.
		plasma_desc_destroy(&dA);

	//	d_print_mat(m, m, A, lda);

		double time = tvmin;
		double Gflops = 1.0/3.0*m*m*m*1e-9/time;

		printf("%d\t%d\t%e\t%f\n", mm, nb, time, Gflops);

		free(A);

		}

	plasma_finalize();

	return 0;

	}


