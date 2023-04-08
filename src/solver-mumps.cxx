
#include "solver-mumps.h"

#if USE_MUMPS

#define __STDC_VERSION__ 199901L

extern "C" {
#include <dmumps_c.h>
}

#define ICNTL(I)        icntl[(I)-1]
#define CNTL(I)         cntl[(I)-1]
#define INFO(I)         info[(I)-1]
#define INFOG(I)        infog[(I)-1]

typedef struct MPS_DATA_
{
    DMUMPS_STRUC_C     *did;
    bool               factored;
    int                ordering;
    int                analysis_type;
    int                refine;
    int                fillin;
    int                verb;

} MPS_DATA;

static const int ordering_values[] = {7, 0, 2, 3, 4, 5, 6};
static int ordering = 0;
static int analysis_type = 0;
static int refine = 0;
static int fillin = 50;

static int mumps_init(MPS_DATA *mpsd, int verb)
{
    mpsd->ordering             = ordering;
    mpsd->analysis_type        = analysis_type;
    mpsd->refine               = refine;
    mpsd->fillin               = fillin;
    mpsd->verb                 = verb;

    mpsd->did = lssp_calloc<DMUMPS_STRUC_C>(1);
    mpsd->did->par = 1;
    mpsd->did->sym = 0;
    mpsd->did->job = -1;

    dmumps_c(mpsd->did);
    mpsd->factored = false;

    return 0;
}

static int mumps_destroy(MPS_DATA *mpsd)
{
    if (mpsd == NULL) return 0;

    lssp_free(mpsd->did->irn);
    lssp_free(mpsd->did->jcn);
    lssp_free(mpsd->did->a);
    lssp_free(mpsd->did->rhs);
    mpsd->did->job = -2;

    dmumps_c(mpsd->did);

    lssp_free(mpsd->did);
    lssp_free(mpsd);

    return 0;
}

static int mumps_fac(MPS_DATA *mpsd, LSSP_SOLVER *solver)
{
    int i, j, k, nz, *icntl, *job;
    int nrow = solver->A.num_rows;
    double time0;
    lssp_mat_csr A = solver->A;

    if (mpsd->factored) return 0;

    nz = A.Ap[nrow];
    mpsd->did->nz = nz;
    mpsd->did->irn = lssp_malloc<int>(nz);
    mpsd->did->jcn = lssp_malloc<int>(nz);
    mpsd->did->a = lssp_malloc<double>(nz);

    mpsd->did->n = nrow;
    mpsd->did->nrhs = 1;
    mpsd->did->lrhs = nrow;
    mpsd->did->rhs = lssp_malloc<double>(mpsd->did->n);

    /* matrix */
    nz = 0;
    for (i = 0; i < nrow; i++) {
        int end = A.Ap[i + 1];

        for (j = A.Ap[i]; j < end; j++) {
            k = A.Aj[j];
            mpsd->did->a[nz] = A.Ax[j];
            mpsd->did->irn[nz] = 1 + i;
            mpsd->did->jcn[nz++] = 1 + k;
        }
    }

    mpsd->factored = true;
    if (mpsd->verb > 0) lssp_printf("mumps: analyse and factor matrix.\n");

    /* parameters */
    icntl = &mpsd->did->ICNTL(0);
    job = &mpsd->did->job;

    /* no outputs */
    if (solver->verb <= 0) {
        icntl[1] = -1;
        icntl[2] = -1;
        icntl[3] = -1;
    }

    icntl[4] = solver->verb;
    icntl[18] = 0;

    /* ordering */
    icntl[7] = ordering_values[mpsd->ordering];
    icntl[13] = 1;
    icntl[28] = mpsd->analysis_type;

    /* analyse */
    time0 = lssp_get_time();
    *job = 1;
    dmumps_c(mpsd->did);

    if (mpsd->verb > 0) lssp_printf("mumps: analysis: %lg s\n", lssp_get_time() - time0);
    if (mpsd->fillin >= 0) icntl[14] = mpsd->fillin;

    /* factorize */
    mpsd->did->CNTL(2) = 0.;
    icntl[10] = mpsd->refine;
    icntl[21] = 0;

    time0 = lssp_get_time();
    *job = 2;
    dmumps_c(mpsd->did);

    if (mpsd->verb > 0) {
        lssp_printf("mumps: factorize: %lg s\n", lssp_get_time() - time0);
    }

    /* check return code */
    i = mpsd->did->INFOG(1);

    if (i < 0) {
        lssp_printf("mumps: return code of MUMPS is %d\n", i);
        lssp_error(1, "mumps: abort.\n");
    }
    else if (i > 0) {
        lssp_warning("mumps: return a non-zero code: %d\n", i);
    }

    return 0;
}

void lssp_solver_mumps_create(LSSP_SOLVER &s)
{
    /* init */
    s.mumps = lssp_calloc<MPS_DATA>(1);
    mumps_init(s.mumps, s.verb);
}

void lssp_solver_mumps_destroy(LSSP_SOLVER &s)
{
    MPS_DATA *mpsd = s.mumps;

    if (mpsd == NULL) return;

    /* destroy */
    mumps_destroy(mpsd);
    s.mumps = NULL;
}

int lssp_solver_mumps(LSSP_SOLVER *solver)
{
    int i, nrow = solver->A.num_rows;
    double time0;
    double *x;
    MPS_DATA *mpsd = solver->mumps;

    /* assemble */
    mumps_fac(mpsd, solver);

    /* rhs */
    for (i = 0; i < nrow; i++) {
        mpsd->did->rhs[i] = solver->rhs.d[i];
    }

    /* solve */
    time0 = lssp_get_time();
    mpsd->did->job = 3;
    dmumps_c(mpsd->did);

    if (mpsd->verb > 0) {
        lssp_printf("mumps: solve: %lg s\n", lssp_get_time() - time0);
    }

    /* solution */
    x = solver->x.d;
    for (i = 0; i < nrow; i++) x[i] = mpsd->did->rhs[i];

    solver->residual = 0.0;
    return solver->nits = 1;
}

#endif
