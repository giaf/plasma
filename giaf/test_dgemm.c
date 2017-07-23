#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "plasma.h"

int main()
	{

	printf("\nciao init\n");

	int ii;

	// problem description

	plasma_enum_t transa = PlasmaNoTrans;
	plasma_enum_t transb = PlasmaNoTrans;
	char cn = 'n';
	char ct = 't';

	int m = 12;
	int n = 12;
	int k = 12;

	double alpha = 1.0;
	double beta = 0.0;

	int Am, An;
	int Bm, Bn;
	int Cm, Cn;

	if(transa == PlasmaNoTrans)
		{
		Am = m;
		An = k;
		}
	else
		{
		Am = k;
		An = m;
		}
	if(transb == PlasmaNoTrans)
		{
		Bm = k;
		Bn = n;
		}
	else
		{
		Bm = n;
		Bn = k;
		}
	Cm = m;
	Cn = n;

	int lda = Am;
	int ldb = Bm;
	int ldc = Cm;

	double *A = (double*)malloc((size_t)lda*An*sizeof(double));
	for(ii=0; ii<m*k; ii++) 
		A[ii] = ii;

	double *B = (double*)malloc((size_t)ldb*Bn*sizeof(double));
	for(ii=0; ii<k*n; ii++)
		B[ii] = 0.0;
	int Bmin = k<n ? k : n;
	for(ii=0; ii<Bmin; ii++)
		B[ii*(ldb+1)] = 1.0;

	double *C = (double*)malloc((size_t)ldc*Cn*sizeof(double));

	omp_set_num_threads(1);

	openblas_set_num_threads(1);

	plasma_init();

	// plasma tuning parameters
	plasma_set(PlasmaTuning, PlasmaDisabled);
	plasma_set(PlasmaNb, 12);

	int rep;
	int nrep = 1;
	struct timeval tv0, tv1;

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		plasma_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
//		dgemm_(&cn, &cn, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
		}

	gettimeofday(&tv1, NULL); // stop

//	d_print_mat(m, n, C, ldc);

	double time = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);
	double Gflops = 2.0*m*n*k*1e-9/time;

	printf("\ntime   = %e\n", time);
	printf("\nGflops = %f\n", Gflops);

	plasma_finalize();

	free(A);
	free(B);
	free(C);

	printf("\nciao end\n");

	return 0;

	}

