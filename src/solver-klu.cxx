
#include "solver-klu.h"

#if USE_SSPARSE

#include <klu.h>

int lssp_solver_klu(LSSP_SOLVER *solver)
{
    double time0;
    double *x = solver->x.d;
    lssp_mat_csr A = solver->A, Ta;

    klu_symbolic *sym ;
    klu_numeric *num ;
    klu_common Common ;
    int i ;

    /* solve */
    time0 = lssp_get_time();
    Ta = lssp_mat_transpose(A);
    for (i = 0 ; i < A.num_rows; i++) x[i] = solver->rhs.d[i];

    /* default pars */
    klu_defaults(&Common);

    sym = klu_analyze(Ta.num_rows, Ta.Ap, Ta.Aj, &Common);
    num = klu_factor(Ta.Ap, Ta.Aj, Ta.Ax, sym, &Common) ;
    klu_solve(sym, num, Ta.num_rows, 1, x, &Common);

    if (solver->verb > 0) {
        lssp_printf("ssparse, klu: solve: %lg s\n", lssp_get_time() - time0);
    }

    klu_free_symbolic(&sym, &Common);
    klu_free_numeric(&num, &Common) ;
    lssp_mat_destroy(Ta);

    solver->residual = 0.0;
    return solver->nits = 1;
}

#endif
