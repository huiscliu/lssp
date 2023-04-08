
#include "solver-qrmumps.h"

#if USE_QR_MUMPS

extern "C" {
#include <dqrm_mumps.h>
}

int lssp_solver_qr_mumps(LSSP_SOLVER *solver)
{
    int i, nrow = solver->A.num_rows;
    double time0;
    struct dqrm_spmat_type_c A, rA;
    lssp_mat_csr Ai = solver->A;
    double *b, *x;

    /* assemble */
    dqrm_spmat_init_c(&A);

    A.m = A.n = nrow;
    A.nz = solver->A.num_nnzs;
    A.icntl[1] = 0;

    A.irn = (int *)malloc(A.nz * sizeof(int));
    A.jcn = (int *)malloc(A.nz * sizeof(int));
    A.val = (double *)malloc(A.nz * sizeof(double));
    x = (double *)malloc(nrow * sizeof(double));
    b = (double *)malloc(nrow * sizeof(double));

    for (i = 0; i < nrow; i++) {
        int end = Ai.Ap[i + 1];

        for (int j = Ai.Ap[i]; j < end; j++) {
            A.irn[j] = i + 1;
            A.jcn[j] = Ai.Aj[j] + 1;
            A.val[j] = Ai.Ax[j];
        }

        x[i] = solver->x.d[i];
        b[i] = solver->rhs.d[i];
    }

    /* solve */
    dqrm_pseti_c(&A, "qrm_keeph", qrm_yes_);
    dqrm_pseti_c(&A, "qrm_ordering", qrm_auto);

    time0 = lssp_get_time();
    dqrm_analyse_c(&A, 'n');
    if (solver->verb > 0) {
        lssp_printf("qr-mumps: analyse: %lg s\n", lssp_get_time() - time0);
    }

    time0 = lssp_get_time();
    dqrm_factorize_c(&A, 'n');
    if (solver->verb > 0) {
        lssp_printf("qr-mumps: factorize: %lg s\n", lssp_get_time() - time0);
    }

    /* rmat */
    dqrm_spmat_init_c(&rA);
    dqrm_get_r_c(&A, &rA);

    time0 = lssp_get_time();
    dqrm_apply_c(&A, 't', b, 1);
    dqrm_solve_c(&A, 'n', b, x, 1);
    if (solver->verb > 0) {
        lssp_printf("qr-mumps: solve: %lg s\n", lssp_get_time() - time0);
    }

    /* solution */
    for (i = 0; i < nrow; i++) solver->x.d[i] = x[i];

    /* destroy */
    dqrm_spmat_destroy_c(&A);
    free(A.irn);
    free(A.jcn);
    free(A.val);
    free(x);
    free(b);

    solver->residual = 0.0;
    return solver->nits = 1;
}

#endif
