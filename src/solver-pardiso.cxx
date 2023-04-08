
#include "solver-pardiso.h"

#if USE_PARDISO

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

int lssp_solver_pardiso(LSSP_SOLVER * solver)
{
    lssp_mat_csr A = solver->A;

    MKL_INT n = A.num_cols;
    MKL_INT nnz = A.num_nonzeros;
    double *a = A.Ax;
    int i;

    MKL_INT mtype = 11;         /* real unsymmetric matrix */
    MKL_INT nrhs = 1;
    MKL_INT idum;               /* integer dummy. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT *ia;
    MKL_INT *ja;

    double *f = solver->rhs.d;
    double *x = solver->x.d;
    void *pt[64];
    double ddum;
    double time = lssp_get_time();

    /* to Fortran */
    ia = lssp_malloc<MKL_INT>(n + 1);
    ja = lssp_malloc<MKL_INT>(nnz);
    for ( i = 0; i <= n; i++ ) ia[i] = A.Ap[i] + 1;
    for ( i = 0; i < nnz; i++ ) ja[i] = A.Aj[i] + 1;

    /* Initialize. */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
        pt[i] = NULL;
    }

    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */

    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 0;           /* Do not print statistical information in file */
    error = 0;            /* Initialize error flag */

    /* reordering and symbolic factorization. */
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        lssp_printf("lssp: pardiso, symbolic factorization failed %d!\n",
                error);
        exit(1);
    }

    /* numerical factorization. */
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        lssp_printf("lssp: pardiso, numerical factorization failed %d!\n", error);
        exit(1);
    }

    /* back substitution and iterative refinement */
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, 
            &nrhs, iparm, &msglvl, f, x, &error);

    if (error != 0) {
        lssp_printf("lssp: pardiso, solution failed %d!\n", error);
        exit(1);
    }

    /* clean */
    phase = -1;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    lssp_free(ia);
    lssp_free(ja);

    /* stat */
    time = lssp_get_time() - time;
    if (solver->verb >= 2) {
        lssp_printf("lssp: pardiso, total time: %g\n", time);
    }

    solver->residual = 0.;
    return 1;
}

#endif
