#include "plasma.h"
#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#ifdef HAVE_BLASFEO_API
#include "blasfeo_d_aux.h"
#endif

void plasma_dprint_blasfeo(int m, int n, plasma_desc_t A)
{

#ifdef HAVE_BLASFEO_API
	printf("\nprint A\n");

	int ii, jj;
	int m0, n0;

	struct blasfeo_dmat sA;

	for(ii=0; ii<A.mt; ii++)
	{
		m0 = plasma_tile_mmain(A, ii);
		for(jj=0; jj<A.nt; jj++)
		{
			n0 = plasma_tile_nmain(A, jj);
			printf("\ntile index %d %d size %d %d\n", ii, jj, m0, n0);
			blasfeo_create_dmat(m0, n0, &sA, plasma_tile_addr(A, ii, jj));
			blasfeo_print_dmat(m0, n0, &sA, 0, 0);
		}
	}
#endif

	return;

}
